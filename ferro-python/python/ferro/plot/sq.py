"""Structure factor from ``ferro traj sq``: **one row of two panels per csv**.

Left $S^{N}(Q)$ (neutron), right $S^{X}(Q)$ (X-ray). They are quantities of different probes
on different scales, so they get separate panels rather than one. Several csv give N rows × 2.

Only the two totals are drawn. The weighted partials (three columns per pair) are a diagnostic
decomposition that adds up to the total, not drawn here: ``ferro traj sq`` always writes them
into the same csv, select the columns in pandas to look at them.
"""

import matplotlib.pyplot as plt

import ferro
from ferro.plot import common

# ── Settings: the quantities one changes ─────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "wspace": 0.35,
    "hspace": 0.45,

    "lw": 1.0,

    "xlabel": r"$Q$ (\AA$^{-1}$)",
    # Left panel first: the order of the list is the drawing order
    "panels": [
        ("total_neutron", r"$S^{N}(Q)$"),
        ("total_xrd",     r"$S^{X}(Q)$"),
    ],
    "legend_loc": "upper right",
    "legend_ncol": 1,

    "xlim": (0, None),
    "ylim": (None, None),
    "hline": 1.0,           # S(Q) tends to 1 at large Q; None draws no reference line
    "hline_kw": {"color": "0.6", "lw": 0.5, "dashes": (2, 2), "zorder": 0},
}


def build_row(df, axes, ctx, path):
    """Draw one product on the two panels of a row."""
    for ax, (col, ylabel) in zip(axes, CFG["panels"]):
        if col not in df.columns:
            raise ferro.FerroError(
                f"{path} has no '{col}' column: that run's --weighting left it out; "
                "run it again with --weighting both")
        if CFG["hline"] is not None:
            ax.axhline(CFG["hline"], **CFG["hline_kw"])
        for name in ctx.shown(df):
            sub = df[df["file"] == name].sort_values("q")
            ax.plot(sub["q"], sub[col], color=ctx.colors[name], lw=CFG["lw"],
                    label=ctx.label(name))
            ctx.record_curve(path.stem, name, col, sub["q"], sub[col])
        ax.set_xlabel(CFG["xlabel"])
        ax.set_ylabel(ylabel)
        ax.set_xlim(*CFG["xlim"])
        ax.set_ylim(*CFG["ylim"])


def build_figure(frames, ctx):
    nrows = len(frames)
    ncols = len(CFG["panels"])
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(CFG["panel_w"] * ncols, CFG["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=CFG["wspace"], hspace=CFG["hspace"])

    for row, (path, df) in zip(axes, frames):
        build_row(df, row, ctx, path)

    handles, labels = axes[0][0].get_legend_handles_labels()
    axes[0][0].legend(handles, labels, loc=CFG["legend_loc"], ncol=CFG["legend_ncol"])
    return fig


def plot(inputs, outdir=None, suffix=None):
    """``sq[_<suffix>].pdf`` / ``.png`` / ``_data.csv`` in ``outdir``; returns their paths."""
    with common.style():
        frames = common.read_all(inputs)
        ctx = common.context(frames, outdir)
        return common.save(build_figure(frames, ctx), common.product_name("sq", suffix), outdir, ctx)
