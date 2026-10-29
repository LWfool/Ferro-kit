"""Bond-angle distribution from ``ferro traj angle``: one panel per csv, one curve per trajectory.

**One axis, not two as for gr**: the product's ``count`` and ``p`` are two scales of the same
curve (``p`` is ``count`` normalised), not two quantities; on two axes they would be two
identical lines. ``count`` is for comparing bin by bin with dump2analysis
(``scripts/compare_angle.py``), not for a figure.
"""

import ferro
from ferro.plot import common

# ── Settings: the quantities one changes ─────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = one row
    "wspace": 0.35,
    "hspace": 0.45,

    "lw": 1.0,
    "smooth": 0,            # > 1: running mean of P(θ) over that many bins (useful for fine bins)

    "xlabel": r"$\theta$ (deg)",
    "ylabel": r"$P(\theta)$",
    "title_fmt": "{a}--{b}--{c}",
    "legend_loc": "upper left",
    "legend_ncol": 1,

    "xlim": (0, 180),
    "ylim": (0, None),
}


def triplet_of(df):
    """The panel title comes from the data columns. More than one triplet is refused: one panel
    cannot say which angle it shows."""
    cols = ["end_a", "center", "end_c"]
    tri = df[cols].drop_duplicates()
    if len(tri) != 1:
        raise ferro.FerroError(
            f"the product holds {len(tri)} triplets, more than one panel can show: run "
            "`ferro traj angle -a END_A -b CENTRE -c END_C` again")
    row = tri.iloc[0]
    return str(row["end_a"]), str(row["center"]), str(row["end_c"])


def build_panel(path, df, ax, ctx):
    a, b, c = triplet_of(df)
    for name in ctx.shown(df):
        sub = df[df["file"] == name].sort_values("angle")
        y = sub["p"]
        if CFG["smooth"] > 1:
            y = y.rolling(CFG["smooth"], center=True, min_periods=1).mean()
        ax.plot(sub["angle"], y, color=ctx.colors[name], lw=CFG["lw"],
                label=ctx.label(name))
        ctx.record_curve(path.stem, name, "p", sub["angle"], y)

    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel"])
    ax.set_title(CFG["title_fmt"].format(a=a, b=b, c=c))
    ax.set_xlim(*CFG["xlim"])
    ax.set_ylim(*CFG["ylim"])


def build_figure(frames, ctx):
    fig, flat = common.panel_grid(len(frames), CFG)
    for ax, (path, df) in zip(flat, frames):
        build_panel(path, df, ax, ctx)

    common.shared_legend(flat, flat[0], CFG)
    return fig


def plot(inputs, outdir=None, suffix=None):
    """``angle[_<suffix>].pdf`` / ``.png`` / ``_data.csv`` in ``outdir``; returns their paths."""
    with common.style():
        frames = common.read_all(inputs)
        ctx = common.context(frames, outdir)
        return common.save(build_figure(frames, ctx), common.product_name("angle", suffix), outdir, ctx)
