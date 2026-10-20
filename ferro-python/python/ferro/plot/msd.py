r"""MSD(t) from ``ferro traj msd``: one panel per csv, one curve per trajectory.

``msd_all.csv msd_Zn.csv`` → two panels. A product made with ``--fit-range`` adds a dashed line
of the same colour: the fitted line, **drawn over the fit window only**, with $D \pm$ its error
(cm²/s) in the legend.

Nothing is fitted here. Slope, intercept, D and its error are the values ferro writes into the
``[inputs]`` list: this module only draws, the algorithm exists once, in ferro. The shared part
of the header holds only what is common to the batch (window fractions, element selection);
every per-file number is in ``[inputs]``.

``loglog`` is the diagnostic MDAnalysis recommends: on log-log axes the diffusive regime has
slope 1, and the fit window should lie on it. The guide line is anchored at the start of the
first trajectory's fit window.
"""

import math

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

from ferro.plot import common

# ── Settings: the quantities one changes ─────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = one row; an integer wraps
    "wspace": 0.35,
    "hspace": 0.7,          # the legend hangs below each panel and needs the room

    # ferro's time is always fs; only the displayed unit changes here
    "time_unit": "ps",      # "ps" or "fs"
    "lw": 1.0,
    "fit_lw": 1.2,
    "fit_dash": (4, 2),
    "comp_lw": 0.7,
    "comp_styles": {"msd_x": ":", "msd_y": "-.", "msd_z": (0, (1, 3))},
    "ref_color": "0.5",

    "ylabel": r"MSD (\AA$^2$)",
    # The legend hangs below the panel: a label carrying D ± error is wider than a 3.3-inch
    # panel and would cover the curves in any corner inside it (a diffusive system fills the
    # diagonal, a glass the top edge). Each panel has its own D, so it cannot be one figure
    # legend either
    "legend_anchor": (0.5, -0.22),
    "legend_fontsize": 6,

    "xlim": (0, None),      # linear axes only
    "ylim": (0, None),
}

TIME_SCALE = {"fs": 1.0, "ps": 1e-3}


def fit_of(inputs, name):
    """One trajectory's fit; None without ``--fit-range`` or when its row failed."""
    if "slope" not in inputs.columns or name not in inputs.index:
        return None
    row = inputs.loc[name]
    keys = ("t_lo", "t_hi", "slope", "intercept", "d_ang2_per_fs", "d_err", "r2")
    vals = {k: float(row[k]) for k in keys}
    return None if math.isnan(vals["slope"]) else vals


def d_label(fit):
    r"""``$D = (1.75 \pm 1.45)\times10^{-8}$ cm$^2$/s``, the exponent taken from D itself."""
    d = fit["d_ang2_per_fs"] * 0.1          # Å²/fs → cm²/s
    err = fit["d_err"] * 0.1
    exp = math.floor(math.log10(abs(d))) if d != 0 else 0
    m = d / 10 ** exp
    body = f"{m:.2f}" if math.isnan(err) else f"({m:.2f} \\pm {err / 10 ** exp:.2f})"
    return f"$D = {body}\\times10^{{{exp}}}$ cm$^2$/s"


def build_panel(path, df, ax, ctx, components, loglog):
    shared, inputs = common.read_header(path)
    scale = TIME_SCALE[CFG["time_unit"]]
    first_fit = None

    for name in ctx.shown(df):
        sub = df[df["file"] == name].sort_values("time")
        if loglog:
            # MSD = 0 at t = 0, which is -inf on log axes
            sub = sub[sub["time"] > 0]
        c = ctx.colors[name]
        fit = fit_of(inputs, name)
        label = ctx.label(name) if fit is None else f"{ctx.label(name)}, {d_label(fit)}"
        ax.plot(sub["time"] * scale, sub["msd"], color=c, lw=CFG["lw"], label=label)
        ctx.record_curve(path.stem, name, "msd", sub["time"] * scale, sub["msd"])

        if components:
            for col, ls in CFG["comp_styles"].items():
                ax.plot(sub["time"] * scale, sub[col], color=c,
                        lw=CFG["comp_lw"], ls=ls, alpha=0.8)
                ctx.record_curve(path.stem, name, col, sub["time"] * scale, sub[col])

        if fit is not None:
            t = np.linspace(fit["t_lo"], fit["t_hi"], 50)
            # Above the noisy curve, or the dashed line of the same colour drowns in it
            ax.plot(t * scale, fit["slope"] * t + fit["intercept"], color=c,
                    lw=CFG["fit_lw"], dashes=CFG["fit_dash"], zorder=3)
            ctx.record_curve(path.stem, name, "fit", t * scale,
                             fit["slope"] * t + fit["intercept"])
            if first_fit is None:
                first_fit = (fit, sub)

    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
        add_slope_one(ax, first_fit, df, scale)
    else:
        ax.set_xlim(*CFG["xlim"])
        ax.set_ylim(*CFG["ylim"])

    sel = shared.get("elements", "all")
    ax.set_title("all atoms" if sel == "all" else sel)
    ax.set_xlabel(f"$t$ ({CFG['time_unit']})")
    ax.set_ylabel(CFG["ylabel"])
    handles, _ = ax.get_legend_handles_labels()
    if components:
        # The components' line styles join the same legend (a second legend would replace the
        # first): colour is still the trajectory, the grey line style only names the direction
        handles += [mlines.Line2D([], [], color="0.3", lw=CFG["comp_lw"], ls=ls,
                                  label=f"{col.split('_')[1]} component")
                    for col, ls in CFG["comp_styles"].items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=CFG["legend_anchor"],
              fontsize=CFG["legend_fontsize"], frameon=False)


def add_slope_one(ax, first_fit, df, scale):
    """A slope-1 guide spanning at most one decade, never beyond the data's time range.

    Anchor: the start of the first trajectory's fit window — if the window was chosen in the
    diffusive regime, the curve should follow this line. Without a fit, the geometric middle
    of the first trajectory's time axis. An anchor less than a factor 2 from the end moves
    back, or the guide is too short to show a slope.
    """
    if first_fit is not None:
        fit, sub = first_fit
        t0 = fit["t_lo"]
    else:
        sub = df[(df["file"] == common.files_in(df)[0]) & (df["time"] > 0)].sort_values("time")
        t0 = math.sqrt(sub["time"].iloc[0] * sub["time"].iloc[-1])
    t_min, t_max = sub["time"].iloc[0], sub["time"].iloc[-1]
    if min(10 * t0, t_max) < 2 * t0:
        t0 = max(t_max / 10, t_min)
    t1 = min(10 * t0, t_max)
    y0 = float(np.interp(t0, sub["time"], sub["msd"]))
    if not (t0 > 0 and y0 > 0 and t1 > t0):
        return
    ax.plot(np.array([t0, t1]) * scale, [y0, y0 * t1 / t0],
            color=CFG["ref_color"], lw=0.8, ls="--")
    ax.annotate("slope 1", (t1 * scale, y0 * t1 / t0), fontsize=6,
                color=CFG["ref_color"], ha="right", va="bottom")


def build_figure(frames, ctx, components, loglog):
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
        build_panel(path, df, ax, ctx, components, loglog)
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig


def plot(inputs, outdir=None, suffix=None, components=False, loglog=False):
    """``msd[_<suffix>].pdf`` / ``.png`` / ``_data.csv`` in ``outdir``; returns their paths.

    ``components`` adds the x / y / z components (they sum to the total); ``loglog`` draws
    log-log axes with a slope-1 guide, to check that the fit window is in the diffusive regime.
    """
    with common.style():
        frames = common.read_all(inputs)
        ctx = common.context(frames, outdir)
        fig = build_figure(frames, ctx, components, loglog)
        return common.save(fig, common.product_name("msd", suffix), outdir, ctx)
