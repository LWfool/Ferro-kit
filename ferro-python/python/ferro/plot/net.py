r"""Network topology from ``ferro net``.

Two kinds of figure of the same shape: **x = composition (the ``file`` column), 100 % stacked
bars, bands = site types**. Both ask one question — how a site's distribution changes with
composition — so they share one template.

    qn      network_qn.csv            one panel per former, bands = Q0…Q4
    cn      network_coordination.csv  one panel per element, bands = coordination number
    bridge  network_linkage.csv       two panels per ligand: the four ligand species
                                      (O_f / O_n / O_b / O_t) and the former–former
                                      linkages (Al-O-Al / Al-O-P / P-O-P…)
    linkmap network_linkage.csv       heat maps, rows = compositions, columns = former pairs
                                      (3 columns for 2 formers, 6 for 3), cells = the states
                                      of the two ends (Al4/5/6 × Q0…Q3)

Stacked bars, unlike lines, draw "sums to 1" as a graphical constraint: the reader sees one
band grow as another shrinks without adding up in their head. The price is that **small
fractions show no trend** — a 2 % band is a thin line. To follow a small fraction, plot that
column alone as a line (or redraw from the exported ``<name>_data.csv``).

**Counting**: ``qn`` counts only **homonuclear** linkages (P–O–P), the $n$ of the literature's
$Q^n_m$. The heteronuclear split (``m_<X>`` columns of ``network_qn_partner.csv``) is not drawn:
two-level stacked bars, hue within lightness, proved hard to read in print, and ``--partner``
was removed (user decision, 2026-09-28). ``ferro net`` still writes that table; filter
$Q^n(m\mathrm{Al})$ in pandas.

## bridge

- **The ligand species** come from the ``network_composition.csv`` beside it (derived from the
  linkage file's name; only the linkage table is given). Without it only the linkage panel
  is drawn
- **The linkages** come from ``network_linkage.csv``, summed over ``(elem_a, elem_b)``.
  ``ferro net`` has already split a triclustered ligand into $C(k,2)$ pair linkages
  (``n_formers`` holds k); they are **counted as such** — the same convention as Qn's
  "by linkage"; their share shows in the ``O_t`` band of the species panel
- The denominator holds only former–former linkages (non-bridging ligands are not in it), and
  the fraction is **summed counts divided**: once rows are merged the mean of per-frame
  fractions cannot be recovered and ``sd`` does not add, so there are no error bars
- With a single former in all inputs the linkage panel is always 100 % and is skipped

## linkmap

- The number on an axis: **a Qn former takes $n$ (homonuclear linkages, P's P–O–P count), any
  other its coordination number** — the same source as the digits of ``ferro net``'s labels.
  Which formers are Qn formers is read from the ``mean_qn`` column of ``[inputs]``, not guessed
  from ``--qn``. Note that Q0 in the Al–P map is "a P bonded to Al with no P–O–P"
- Rows are ``elem_a``, columns ``elem_b`` (linkage's canonical order: Al–P has Al on the rows).
  A map of one element with itself shows only the upper triangle: a bridge has no direction
  and is stored once, and a mirror would make the off-diagonal linkages look doubled
- **Each map is normalised to 100 % on its own** (one decimal in the cells); the colour scale
  is shared down a column, for comparing compositions
- Axis values are the union over all compositions **per map and per end**, so the
  compositions' grids line up; they are not shared across maps by element — Q0 cannot occur in
  the P–O–P map (a P in a P–O–P has n ≥ 1), and sharing would add an empty row. Light grey =
  no composition has it; 0 written = another composition has it, this one does not; a whole
  grey map = this composition has no such linkage
- Triclusters counted as in bridge; nothing is drawn with a single former

## Several csv

The ``-s`` batch tag of ``ferro net`` is the last segment of the name, so
``network_qn_CMD.csv`` and ``network_qn_MLMD.csv`` tell themselves apart. Several csv give
**several bars side by side on one composition axis**, the group named after the suffix
(renamed or hidden in ``labels.csv`` like any key).
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

import ferro
from ferro.plot import common

# ── Settings: the quantities one changes ─────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = one row
    "wspace": 0.35,
    "hspace": 0.55,

    "bar_width": 0.55,      # the width all bars on one composition tick take together
    "bar_gap": 0.08,        # gap between groups side by side (fraction of bar_width)
    "edge_lw": 0.4,         # outline of each band
    "edge_color": "white",

    "xlabel": "",           # the compositions are on the ticks already
    "ylabel": r"Fraction (\%)",
    "title_qn": "{former}",
    "title_cn": "{element}",
    "xtick_rotation": 30,
    # With several csv the group names (CMD / MLMD) stand vertically above the bars; the title
    # has to make room
    "title_pad_grouped": 22,
    "legend_loc": "upper center",
    # The legend hangs below **the whole figure**, not below one panel: the compositions are
    # long rotated labels and would pierce a legend hung on an axis. y leaves room for them
    "legend_bbox": (0.5, -0.16),
    "legend_ncol": 6,
    # bridge's two panels have different categories; each legend hangs to its panel's right
    "panel_legend_bbox": (1.02, 1.0),
    "wspace_panel_legend": 0.9,
    "title_species": "{ligand} species",
    "title_bridge": r"$X$--{ligand}--$X$ linkages",

    "ylim": (0, 100),

    # linkmap
    "map_cell": 0.36,       # side of one cell (inches)
    "map_pad_w": 0.55,      # width / height of a map beyond its cells (ticks and labels)
    "map_pad_h": 0.45,
    "map_wspace": 0.5,      # fraction of the mean panel width
    "map_hspace": 0.3,
    "map_cmap": "Blues",
    "map_fontsize": 5.5,
    "map_text_dark": 0.6,   # white text above this fraction of the column's maximum
    "map_never": "0.88",    # a cell no composition has
    "map_absent": "0.93",   # this composition has no such linkage: the whole map
    "map_title": "{a}--{lig}--{b}",
    "map_cbar_label": r"Fraction (\%)",
}

# Source of the notes that concern the batch rather than one file
PROGRAM = "ferro-py"


# ── Product names ────────────────────────────────────────────────────────────
# Both functions below take ferro's product names apart: `<mode>[_<table>]...[_<-s>]`, built by
# `product_name` in ferro-cli/src/batch.rs. Two readers of one Rust rule, with no binding of it
# in Python; a change of that rule must be followed here (scan E-13).

TABLES = ("network_qn_partner", "network_qn", "network_coordination",
          "network_composition", "network_ligand_type", "network_linkage")


def suffix_of(path):
    """The batch name in a product's file name, the group label of several csv side by side.

    ``network_qn_CMD.csv`` → ``CMD``. Without one (``network_qn.csv``, the usual run without
    ``-s``) the whole stem: a single file needs no tag made up for a figure.
    """
    stem = path.stem
    for table in TABLES:
        if stem.startswith(table + "_"):
            return stem[len(table) + 1:]
    return stem


def composition_of(path):
    """The composition product in the linkage product's directory, same suffix; None if absent."""
    sib = path.with_name(path.name.replace("network_linkage", "network_composition", 1))
    if sib == path or not sib.is_file():
        return None
    return common.read(sib)


# ── Shaping the data ─────────────────────────────────────────────────────────

def series_qn(df, former):
    """[(band label, category value, the band's rows)], by ascending qn (stacked bottom-up)."""
    sub = df[df["former"] == former]
    return [(f"Q{int(q)}", int(q), sub[sub["qn"] == q])
            for q in sorted(sub["qn"].unique())]


def series_cn(df, element):
    """The band label is the coordination number alone, without the element.

    Colours are **assigned by category value across the whole figure**, so there is one shared
    legend and the element is in the panel title. A label ``Al_4`` would hang ``Al`` on the
    same orange in the Zn panel.
    """
    sub = df[df["element"] == element]
    return [(f"CN = {int(c)}", int(c), sub[sub["cn"] == c])
            for c in sorted(sub["cn"].unique())]


# Stacking order of the four ligand species: free → non-bridging → bridging → tricluster
# (increasing number of formers bonded)
SPECIES = ("f", "n", "b", "t")


def by_count(df, key_cols):
    """Sum count over (file, key_cols); fraction = share of the count **within the file**.

    After merging rows this is the only way: the csv's fraction is the mean of per-frame
    fractions, and a sum of rows' means is still a mean, but not once the denominator changes
    (only some rows taken).
    """
    g = df.groupby(["file", *key_cols], sort=False, as_index=False)["count"].sum()
    g["fraction"] = g["count"] / g.groupby("file")["count"].transform("sum")
    return g


def series_species(comp, ligand):
    """The four ligand species, from the composition table's ``<ligand>_<f|n|b|t>`` rows."""
    tags = {f"{ligand}_{t}": t for t in SPECIES}
    sub = comp[(comp["element"] == ligand) & comp["label"].isin(tags)].copy()
    sub["tag"] = sub["label"].map(tags)
    g = by_count(sub, ["tag"])
    return [(rf"{ligand}$_\mathrm{{{t}}}$", i, g[g["tag"] == t])
            for i, t in enumerate(SPECIES) if (g["tag"] == t).any()]


def series_bridge(link, ligand):
    """Former–former linkages summed over (elem_a, elem_b); a tricluster's C(k,2) pairs count."""
    g = by_count(link[link["ligand"] == ligand], ["elem_a", "elem_b"])
    pairs = sorted(set(zip(g["elem_a"], g["elem_b"])))
    return [(f"{a}-{ligand}-{b}", (a, b), g[(g["elem_a"] == a) & (g["elem_b"] == b)])
            for a, b in pairs]


def nontrivial_elements(df):
    """An element with a single coordination number is skipped: a flat 100 % is noise, not
    information.

    The same principle as ``ferro net`` with an empty Qn table: nothing that was measured but
    has no content.
    """
    keep, skipped = [], []
    for e in dict.fromkeys(df["element"]):
        (keep if df[df["element"] == e]["cn"].nunique() > 1 else skipped).append(e)
    return keep, skipped


# ── Drawing ──────────────────────────────────────────────────────────────────

def plain(tex):
    r"""A LaTeX label of the figure → plain text in the export: ``O$_\mathrm{f}$`` → ``O_f``.

    The export is for redrawing; the reader should get names to filter on, not typesetting.
    """
    return (tex.replace("$", "").replace(r"\mathrm{", "").replace("{", "")
               .replace("}", "").replace("--", "-"))


def build_panel(ax, groups, title, colors, handles, seen, ctx):
    """groups: [(group name, segs)], segs from series_*.

    One x tick = one composition; one bar per group on a tick; the bands stacked bottom-up.

    ``colors`` is {category value: colour} **for the whole figure**, and ``handles`` / ``seen``
    accumulate over the whole figure too: a colour must mean the same in every panel, or a
    shared legend lies.
    """
    files = ctx.names
    n_grp = len(groups)
    total_w = CFG["bar_width"]
    w = total_w / n_grp * (1 - CFG["bar_gap"])
    offsets = [(-total_w / 2 + total_w / n_grp * (i + 0.5)) for i in range(n_grp)]

    for gi, (gname, segs) in enumerate(groups):
        xs = [i + offsets[gi] for i in range(len(files))]
        bottom = [0.0] * len(files)
        for label, key, rows in segs:
            vals = [float(rows[rows["file"] == f]["fraction"].sum()) * 100 for f in files]
            ax.bar(xs, vals, width=w, bottom=bottom, color=colors[key],
                   edgecolor=CFG["edge_color"], linewidth=CFG["edge_lw"], zorder=2)
            bottom = [b + v for b, v in zip(bottom, vals)]
            for f, v in zip(files, vals):
                if (rows["file"] == f).any():
                    ctx.record(panel=plain(title), group=gname, file=f, category=plain(label),
                               percent=round(v, 6), count=int(rows[rows["file"] == f]["count"].sum()))
            if key not in seen:
                seen.add(key)
                handles.append((key, Patch(facecolor=colors[key], label=label)))

        if n_grp > 1:
            ax.text(offsets[gi], 101, ctx.label(gname), ha="center", va="bottom",
                    fontsize=plt.rcParams["font.size"] - 1, rotation=90)

    ax.set_xticks(range(len(files)))
    ax.set_xticklabels([ctx.label(f) for f in files],
                       rotation=CFG["xtick_rotation"], ha="right")
    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel"])
    ax.set_ylim(*CFG["ylim"])
    ax.set_title(title, pad=CFG["title_pad_grouped"] if n_grp > 1 else None)
    ax.tick_params(top=False)          # mirrored ticks mean nothing on a categorical axis


def build_figure(panels, groups_of, ctx, panel_legend=False, palette_keys=None):
    """panels: [(title, key)], groups_of(key) -> [(group name, segs)].

    By default one shared legend below the whole figure, colours assigned by category value
    across the figure. ``panel_legend`` is for figures whose panels have different categories
    (bridge): colours per panel, the legend to each panel's right — otherwise the two panels'
    categories together outrun the colour cycle and one colour appears twice in one legend.

    ``palette_keys(key)`` gives **every possible** category of a panel; when given, colours
    follow it rather than the categories present in this batch — otherwise one missing category
    (no O_f) shifts a whole row of colours, and one category differs in colour between figures.
    """
    fig, flat = common.panel_grid(len(panels), CFG,
                                  wspace=CFG["wspace_panel_legend"] if panel_legend else None)

    per_panel = [groups_of(key) for _, key in panels]

    def keys_of(groups_list):
        return sorted({k for groups in groups_list for _, segs in groups for _, k, _ in segs})

    shared = common.color_map(keys_of(per_panel)) if not panel_legend else None
    handles, seen = [], set()
    for ax, (title, key), groups in zip(flat, panels, per_panel):
        if panel_legend:
            handles, seen = [], set()
        fixed = palette_keys(key) if palette_keys else None
        colors = shared or common.color_map(fixed or keys_of([groups]))
        build_panel(ax, groups, title, colors, handles, seen, ctx)
        if panel_legend:
            # Top-down is the reverse of the bottom-up stacking; reversed, the legend follows
            # the bars
            hs = [h for _, h in sorted(handles, key=lambda kh: kh[0], reverse=True)]
            ax.legend(handles=hs, loc="upper left", bbox_to_anchor=CFG["panel_legend_bbox"],
                      frameon=False)
    if not panel_legend:
        # Legend sorted by category value, the order of the bottom-up stacking
        handles = [h for _, h in sorted(handles, key=lambda kh: kh[0])]
        fig.legend(handles=handles, loc=CFG["legend_loc"],
                   bbox_to_anchor=CFG["legend_bbox"], ncol=CFG["legend_ncol"],
                   bbox_transform=fig.transFigure)
    return fig


# ── The four figures ─────────────────────────────────────────────────────────

def run_qn(frames, ctx, elements):
    formers = list(dict.fromkeys(f for _, df in frames for f in df["former"]))

    def groups_of(former):
        return [(suffix_of(p), series_qn(df, former)) for p, df in frames]

    return build_figure([(CFG["title_qn"].format(former=f), f) for f in formers],
                        groups_of, ctx)


def run_cn(frames, ctx, elements):
    if not elements:
        elements, skipped = [], []
        for _, df in frames:
            k, s = nontrivial_elements(df)
            elements += [e for e in k if e not in elements]
            skipped += [e for e in s if e not in skipped]
        for e in skipped:
            if e not in elements:
                common.note(PROGRAM, f"{e} has a single coordination number; skipped "
                            "(it would be a flat 100 %)")
    if not elements:
        raise ferro.FerroError("no element has more than one coordination number; "
                               "name the elements to draw with --element")

    def groups_of(element):
        return [(suffix_of(p), series_cn(df, element)) for p, df in frames]

    return build_figure([(CFG["title_cn"].format(element=e), e) for e in elements],
                        groups_of, ctx)


def run_bridge(frames, ctx, elements):
    comps = {p: composition_of(p) for p, _ in frames}
    for p, c in comps.items():
        if c is None:
            common.note(p, "no network_composition product beside it; the ligand species "
                        "are not drawn")
    ligands = list(dict.fromkeys(l for _, df in frames for l in df["ligand"]))
    formers = {e for _, df in frames for col in ("elem_a", "elem_b") for e in df[col]}

    panels = []
    for lig in ligands:
        if any(c is not None for c in comps.values()):
            panels.append((CFG["title_species"].format(ligand=lig), ("species", lig)))
        if len(formers) > 1:
            panels.append((CFG["title_bridge"].format(ligand=lig), ("bridge", lig)))
    if len(formers) <= 1:
        common.note(PROGRAM, f"a single former in all inputs ({', '.join(sorted(formers))}): "
                    "the linkages are always 100 %, not drawn")
    if not panels:
        raise ferro.FerroError("nothing to draw: no composition product and a single former")

    def groups_of(key):
        what, lig = key
        out = []
        for p, df in frames:
            if what == "species":
                c = comps[p]
                segs = [] if c is None else series_species(c, lig)
            else:
                segs = series_bridge(df, lig)
            out.append((suffix_of(p), segs))
        return out

    def palette_keys(key):
        return list(range(len(SPECIES))) if key[0] == "species" else None

    return build_figure(panels, groups_of, ctx, panel_legend=True, palette_keys=palette_keys)


def qn_formers(frames):
    """The elements in the ``mean_qn`` column of ``[inputs]`` (``P=0.95 Si=1.2``) = the Qn
    formers."""
    out = set()
    for p, _ in frames:
        _, inputs = common.read_header(p)
        if "mean_qn" in inputs.columns:
            for cell in inputs["mean_qn"].dropna():
                out |= {tok.split("=")[0] for tok in str(cell).split() if "=" in tok}
    return out


def link_states(df, qn_set):
    """Add each row's two ends as numbers on the map axes: qn for a Qn former, cn otherwise."""
    df = df.copy()
    for end in ("a", "b"):
        use_qn = df[f"elem_{end}"].isin(qn_set)
        df[f"v_{end}"] = np.where(use_qn, df[f"qn_{end}"], df[f"cn_{end}"]).astype(int)
    same = df["elem_a"] == df["elem_b"]
    # A map of one element with itself shows the upper triangle: the canonical order sorts by
    # (element, qn, cn) while the axis may use cn, so rows with v_a > v_b exist (Al: small qn,
    # large cn) and are folded into the upper triangle
    lo = np.minimum(df["v_a"], df["v_b"])
    hi = np.maximum(df["v_a"], df["v_b"])
    df.loc[same, "v_a"] = lo[same]
    df.loc[same, "v_b"] = hi[same]
    return df


def tick(elem, v, qn_set):
    return f"Q{v}" if elem in qn_set else f"{elem}{v}"


def run_linkmap(frames, ctx, elements):
    qn_set = qn_formers(frames)
    frames = [(p, link_states(df, qn_set)) for p, df in frames]
    formers = {e for _, df in frames for col in ("elem_a", "elem_b") for e in df[col]}
    if len(formers) <= 1:
        common.note(PROGRAM, f"a single former in all inputs ({', '.join(sorted(formers))}): "
                    "no linkmap drawn")
        return None

    # Columns = former pairs, by (ligand, elem_a, elem_b); axis values = union over compositions
    pairs = sorted({(l, a, b) for _, df in frames
                    for l, a, b in zip(df["ligand"], df["elem_a"], df["elem_b"])})
    # axis[(pair, end)]: the values that end of that map takes in any composition
    axis = {}
    for _, df in frames:
        for pair, sel in df.groupby(["ligand", "elem_a", "elem_b"]):
            for end in ("a", "b"):
                axis.setdefault((pair, end), set()).update(int(v) for v in sel[f"v_{end}"])
    axis = {k: sorted(vs) for k, vs in axis.items()}
    if any(a == b for _, a, b in pairs):
        # A map of one element with itself uses one axis for both ends; only then does the
        # upper triangle mean something
        for pair in pairs:
            if pair[1] == pair[2]:
                both = sorted(set(axis[(pair, "a")]) | set(axis[(pair, "b")]))
                axis[(pair, "a")] = axis[(pair, "b")] = both

    # Rows = (group, composition); with several groups the row name carries the group
    rows = [(p, suffix_of(p), f) for p, df in frames for f in ctx.shown(df)]
    multi = len(frames) > 1

    # The count matrix per (row, column)
    mats = {}
    for p, df in frames:
        for f in ctx.shown(df):
            sub = df[df["file"] == f]
            for pair in pairs:
                lig, a, b = pair
                ya, xb = axis[(pair, "a")], axis[(pair, "b")]
                m = np.zeros((len(ya), len(xb)))
                sel = sub[(sub["ligand"] == lig) & (sub["elem_a"] == a) & (sub["elem_b"] == b)]
                for va, vb, c in zip(sel["v_a"], sel["v_b"], sel["count"]):
                    m[ya.index(va), xb.index(vb)] += c
                mats[(p, f, pair)] = m

    nr, nc = len(rows), len(pairs)
    # Width ratios by cell count only: with aspect=equal that makes every map's cells equal
    cols_x = [len(axis[(pair, "b")]) for pair in pairs]
    rows_y = max(len(axis[(pair, "a")]) for pair in pairs)
    fig_w = sum(cols_x) * CFG["map_cell"] + nc * CFG["map_pad_w"] + 0.4
    fig_h = nr * (rows_y * CFG["map_cell"] + CFG["map_pad_h"]) + 0.7
    fig, axes = plt.subplots(nr, nc, figsize=(fig_w, fig_h), squeeze=False,
                             gridspec_kw={"width_ratios": cols_x})
    fig.subplots_adjust(wspace=CFG["map_wspace"], hspace=CFG["map_hspace"])
    cmap = mpl.colormaps[CFG["map_cmap"]]

    for j, pair in enumerate(pairs):
        lig, a, b = pair
        pcts = {}
        for p, g, f in rows:
            m = mats[(p, f, pair)]
            tot = m.sum()
            pcts[(p, f)] = m / tot * 100 if tot else None
        # Cells no composition has → light grey
        seen = sum((mats[(p, f, pair)] > 0) for p, _, f in rows) > 0
        vmax = max((v.max() for v in pcts.values() if v is not None), default=1.0) or 1.0
        norm = mpl.colors.Normalize(0, vmax)
        same = a == b

        for i, (p, g, f) in enumerate(rows):
            ax = axes[i][j]
            pct, cnt = pcts[(p, f)], mats[(p, f, pair)]
            ny, nx = cnt.shape
            rgba = np.ones((ny, nx, 4))
            for y in range(ny):
                for x in range(nx):
                    if same and y > x:
                        continue                      # lower triangle left blank
                    if pct is None:
                        rgba[y, x] = mpl.colors.to_rgba(CFG["map_absent"])
                    elif not seen[y, x]:
                        rgba[y, x] = mpl.colors.to_rgba(CFG["map_never"])
                    else:
                        v = pct[y, x]
                        rgba[y, x] = cmap(norm(v))
                        # Non-zero but 0.0 to one decimal is written <0.1, apart from a real 0
                        txt = "0" if v == 0 else (f"{v:.1f}" if v >= 0.05 else "<0.1")
                        ax.text(x, y, txt, ha="center", va="center",
                                fontsize=CFG["map_fontsize"],
                                color="white" if v > CFG["map_text_dark"] * vmax else "black")
                        ctx.record(panel=f"{a}-{lig}-{b}", group=g, file=f,
                                   row=tick(a, axis[(pair, "a")][y], qn_set),
                                   col=tick(b, axis[(pair, "b")][x], qn_set),
                                   percent=round(float(v), 6), count=int(cnt[y, x]))
            ax.imshow(rgba, aspect="equal")
            if pct is None:
                ax.text((nx - 1) / 2, (ny - 1) / 2, "absent", ha="center", va="center",
                        fontsize=CFG["map_fontsize"] + 1, color="0.4")
            ax.set_xticks(range(nx))
            ax.set_xticklabels([tick(b, v, qn_set) for v in axis[(pair, "b")]])
            ax.set_yticks(range(ny))
            ax.set_yticklabels([tick(a, v, qn_set) for v in axis[(pair, "a")]])
            ax.tick_params(which="both", length=0, top=False, right=False)
            ax.minorticks_off()
            ax.tick_params(labelsize=CFG["map_fontsize"] + 0.5)
            if i == 0:
                ax.set_title(CFG["map_title"].format(a=a, lig=lig, b=b))

        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[-1][j],
                          orientation="horizontal", location="bottom", shrink=0.9, pad=0.12,
                          aspect=15)
        cb.set_label(CFG["map_cbar_label"], fontsize=CFG["map_fontsize"] + 0.5)
        cb.ax.tick_params(labelsize=CFG["map_fontsize"], top=False)
        cb.ax.minorticks_off()
    # A colorbar borrowing room from the last row resets the anchor; align all tops at the end
    for ax in axes.ravel():
        ax.set_anchor("N")
    # Row names are centred on **the whole row's slot**, not hung on the first column's ylabel:
    # the first column may be a small map of two rows, and a longer name would leave the canvas
    left = mpl.transforms.offset_copy(fig.transFigure, fig=fig, x=-30, units="points")
    for i, (p, g, f) in enumerate(rows):
        slot = axes[i][0].get_subplotspec().get_position(fig)
        name = ctx.label(f) if not multi else f"{ctx.label(g)}: {ctx.label(f)}"
        fig.text(slot.x0, slot.y1 - (rows_y * CFG["map_cell"]) / fig_h / 2, name,
                 transform=left, rotation=90, ha="center", va="center")
    return fig


RUN = {"qn": run_qn, "cn": run_cn, "bridge": run_bridge, "linkmap": run_linkmap}

# Columns each figure reads; checked up front so the wrong table names itself
# instead of surfacing as a bare KeyError deep inside a plot
NEEDS = {
    "qn": {"file", "former", "qn", "count"},
    "cn": {"file", "element", "cn", "count"},
    "bridge": {"file", "ligand", "elem_a", "elem_b", "count"},
    "linkmap": {"file", "ligand", "elem_a", "elem_b", "qn_a", "qn_b", "cn_a", "cn_b", "count"},
}


def plot(inputs, table, outdir=None, suffix=None, elements=None):
    """``net_<table>[_<suffix>].pdf`` / ``.png`` / ``_data.csv`` in ``outdir``; returns their
    paths (none when there is nothing to draw).

    ``table`` is ``qn``, ``cn``, ``bridge`` or ``linkmap``; ``elements`` limits ``cn`` to these
    elements (default: those with more than one coordination number).
    """
    if table not in RUN:
        raise ferro.FerroError(f"unknown network figure '{table}'; one of {', '.join(RUN)}")
    with common.style():
        frames = common.read_all(inputs)
        for p, df in frames:
            missing = NEEDS[table] - set(df.columns)
            if missing:
                raise ferro.FerroError(f"{p}: '{table}' needs column(s) {', '.join(sorted(missing))}; "
                                       "is this the right ferro net table?")

        # Group names show only with several csv side by side; with one csv the name falls back
        # to the whole stem and registering it would only clutter labels.csv
        groups = [suffix_of(p) for p, _ in frames] if len(frames) > 1 else []
        ctx = common.context(frames, outdir, groups=groups)
        # A group with show=0 in labels.csv is left out whole
        frames = [(p, df) for p, df in frames if not groups or ctx.group_shown(suffix_of(p))]
        if not frames:
            raise ferro.FerroError(f"{ctx.path} sets show=0 for every group of this run: "
                                   "nothing to draw")

        fig = RUN[table](frames, ctx, elements)
        if fig is None:
            return []
        return common.save(fig, common.product_name(f"net_{table}", suffix), outdir, ctx)
