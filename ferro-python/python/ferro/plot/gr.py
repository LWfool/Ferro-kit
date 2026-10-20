"""g(r) / CN(r) from ``ferro traj gr``: one panel per csv, one curve per trajectory.

``gr_P-O.csv gr_Al-O.csv`` → two panels, each overlaying the compositions. g(r) solid on the
left axis, CN(r) dashed on the right, one colour per trajectory.
"""

import matplotlib.pyplot as plt

import ferro
from ferro.plot import common

# ── Settings: the quantities one changes ─────────────────────────────────────

CFG = {
    # One panel (inches); the figure width is panel_w × columns
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = one row; an integer wraps
    "wspace": 0.55,         # the right axis label needs the room; do not go much lower
    "hspace": 0.45,

    "gr_lw": 1.0,
    "cn_lw": 0.9,
    "cn_dash": (3, 1.5),    # dash pattern of CN

    "xlabel": r"$r$ (\AA)",
    "ylabel_gr": r"$g(r)$",
    "ylabel_cn": r"$CN(r)$",
    "title_fmt": "{center}--{neighbor}",
    "legend_loc": "upper right",
    "legend_ncol": 1,

    "xlim": (0, None),      # None = automatic
    "gr_ylim": (0, None),
    "cn_ylim": (0, None),
}


def pair_of(df):
    """The panel title comes from the **data columns**, not the file name.

    A file may have been renamed; the ``center`` / ``neighbor`` columns are not. A gr product
    usually holds one pair (``-a P -b O``); ``gr_all.csv`` from a run without a selection holds
    every ordered pair, and a title cannot say which one is meant, so it is refused.
    """
    pairs = df[["center", "neighbor"]].drop_duplicates()
    if len(pairs) != 1:
        raise ferro.FerroError(
            f"the product holds {len(pairs)} pairs, more than one panel can show: run "
            "`ferro traj gr -a CENTRE -b NEIGHBOUR` again, or select the rows in pandas")
    row = pairs.iloc[0]
    return str(row["center"]), str(row["neighbor"])


def build_panel(path, df, ax, ctx):
    """Draw one product on ax: g(r) solid on the left axis, CN(r) dashed on the right."""
    center, neighbor = pair_of(df)

    # The science style mirrors the main axis' ticks on the right, where they would collide
    # with the CN axis
    ax.tick_params(right=False)
    ax_cn = ax.twinx()
    ax_cn.tick_params(right=True)

    for name in ctx.shown(df):
        sub = df[df["file"] == name].sort_values("r")
        c = ctx.colors[name]
        ax.plot(sub["r"], sub["gr"], color=c, lw=CFG["gr_lw"], label=ctx.label(name))
        ax_cn.plot(sub["r"], sub["cn"], color=c, lw=CFG["cn_lw"],
                   dashes=CFG["cn_dash"])
        ctx.record_curve(path.stem, name, "gr", sub["r"], sub["gr"])
        ctx.record_curve(path.stem, name, "cn", sub["r"], sub["cn"])

    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel_gr"])
    ax_cn.set_ylabel(CFG["ylabel_cn"])
    ax.set_title(CFG["title_fmt"].format(center=center, neighbor=neighbor))
    ax.set_xlim(*CFG["xlim"])
    ax.set_ylim(*CFG["gr_ylim"])
    ax_cn.set_ylim(*CFG["cn_ylim"])
    return ax_cn


def build_figure(frames, ctx):
    """frames: [(path, df)]; returns the figure."""
    n = len(frames)
    ncols = CFG["ncols"] or n
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(CFG["panel_w"] * ncols, CFG["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=CFG["wspace"], hspace=CFG["hspace"])

    flat = axes.ravel()
    for ax, (path, df) in zip(flat, frames):
        build_panel(path, df, ax, ctx)
    for ax in flat[n:]:
        ax.set_visible(False)

    # One legend: the colours mean the same in every panel
    handles, labels = flat[0].get_legend_handles_labels()
    flat[0].legend(handles, labels, loc=CFG["legend_loc"], ncol=CFG["legend_ncol"])
    return fig


def plot(inputs, outdir=None, suffix=None):
    """``gr[_<suffix>].pdf`` / ``.png`` / ``_data.csv`` in ``outdir``; returns their paths."""
    with common.style():
        frames = common.read_all(inputs)
        ctx = common.context(frames, outdir)
        return common.save(build_figure(frames, ctx), common.product_name("gr", suffix), outdir, ctx)
