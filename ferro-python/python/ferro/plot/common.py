"""The shared layer of the figure modules: style, reading ferro's csv, ``labels.csv``, export.

The comparison scripts in ``scripts/`` (``ferrocmp.py``) also read ferro's csv, but they
compare against reference programs and draw only to show differences; after reading the two
share no code, so they stay apart.

## Style

``['science', 'vibrant']`` + LaTeX. ``science`` is a four-sided frame with inward ticks,
``vibrant`` Paul Tol's colour-blind-safe sequence of 7. Two findings from testing:

- **Strings need no escaping.** matplotlib's usetex preamble loads the ``underscore`` package
  (``[strings]`` option), so ``43Z43P15A_NPT_5``, ``P_3``, ``O_b`` and ``50% Zn`` all pass.
- The only characters that still make **the whole figure fail to render** are ``&``, ``#`` and
  ``^``. Element symbols and site labels cannot hold them (ferro's file names allow only
  ``[A-Za-z0-9_+-]``); only a custom file name can. Rename that file rather than adding an
  escaping layer here: every string going into a figure would have to pass through it, and one
  missed is the same failure, only harder to find.

## Colours

Colours are **always assigned by the ``file`` column** (one trajectory, one colour), the same
across panels and figures. That is what all figures share: ``file`` is the composition or
system, the dimension the reader follows.

## Display names: ``labels.csv``

One ``labels.csv`` in the output directory, **shared by every figure** (columns
``key,label,show``). A composition's name is the same throughout a paper and should not be
changed once in the g(r) figure and again in the Qn figure.

- The first run writes it, ``label`` = the key, ``show`` = 1
- **Row order is display order** (x axis, legend, colour assignment); ``show=0`` is not drawn
- On every later run existing rows are left untouched; keys new to this run are appended and
  reported as a FerroNote (a mistyped key shows up as such a new row: that report is the guard
  instead of an error); keys in the file but not in this run's data are ignored silently —
  they belong to another batch, not a mistake
- A key is a composition (the ``file`` column) or, with several csv side by side, a group (the
  suffix of the file name)

## Data export: ``<name>_data.csv``

Every figure writes a long table next to the pdf with exactly the numbers drawn (after
smoothing and filtering), to redraw it in a journal's style. Both the original ``file`` and the
mapped ``label`` are columns; row order is drawing order.
"""

import contextlib
import csv
import warnings
from pathlib import Path

import ferro

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scienceplots  # noqa: F401  (importing registers the styles)
except ImportError as e:
    raise ferro.FerroError(f"plotting needs matplotlib, pandas and scienceplots ({e}): "
                           "pip install 'ferro[plot]'") from None

# ── Export settings ──────────────────────────────────────────────────────────

STYLE = ["science", "vibrant"]

#: The png is only for a look; the pdf is the product
PNG_DPI = 600

#: 42 = TrueType, so the text in the pdf can be selected and edited in Illustrator
PDF_FONTTYPE = 42


@contextlib.contextmanager
def style():
    """The plotting style for one figure, without changing the caller's matplotlib settings."""
    with plt.style.context(STYLE), mpl.rc_context({"pdf.fonttype": PDF_FONTTYPE,
                                                   "ps.fonttype": PDF_FONTTYPE}):
        yield


def note(source, detail):
    """A FerroNote in the layout of ferro's other notices (``source: FerroNote: detail``)."""
    w = ferro.FerroNote(f"{source}: {detail}")
    w.source, w.line, w.detail = str(source), None, detail
    warnings.warn(w, stacklevel=2)


def product_name(stem, suffix):
    """``<stem>[_<suffix>]``: the ``-s`` batch tag is the last segment, as in ferro's products."""
    return f"{stem}_{suffix}" if suffix else stem


# ── Reading ──────────────────────────────────────────────────────────────────

def read(path):
    """Read one ferro product.

    The ``#`` block is for people (version, shared parameters, the ``[inputs]`` list), so
    ``comment="#"`` drops it whole. A missing value is an empty field → NaN: **do not count it
    as 0**.
    """
    path = Path(path)
    df = pd.read_csv(path, comment="#")
    if df.empty:
        raise ferro.FerroError(f"{path} has no data rows (only a header?)")
    return df


def read_all(paths):
    """[(path, DataFrame)] in the order given."""
    return [(Path(p), read(p)) for p in paths]


def read_header(path):
    """Read the ``#`` block: the shared parameters (``key = value``) and the ``[inputs]`` list.

    ``[inputs]`` is a fixed-width table for people; ``composition`` / ``species`` and the failure
    reasons contain spaces, so it is cut at the start columns of the header's fields, not at
    whitespace. ``-`` means "no value" → NaN.
    """
    lines = [l[2:].rstrip("\n") for l in open(path, encoding="utf-8") if l.startswith("# ")]
    shared = {}
    inputs = None
    for i, line in enumerate(lines):
        if line == "[inputs]":
            head = lines[i + 1]
            starts = [j for j, ch in enumerate(head)
                      if ch != " " and (j == 0 or head[j - 1] == " ")]
            names = head.split()
            rows = []
            for row in lines[i + 2:]:
                if row.startswith("---"):
                    break
                bounds = starts[1:] + [None]
                rows.append([row[a:b].strip() for a, b in zip(starts, bounds)])
            inputs = pd.DataFrame(rows, columns=names).replace("-", np.nan)
            break
        if "=" in line:
            key, _, val = line.partition("=")
            shared[key.strip()] = val.strip()
    if inputs is None:
        raise ferro.FerroError(f"{path} has no [inputs] list: not a ferro analysis product?")
    return shared, inputs.set_index("file")


def files_in(df):
    """The values of the ``file`` column, **in order of appearance**.

    The order is not alphabetical: it is the order of ``ferro -i`` (``expand_inputs`` keeps it,
    ``concat_union`` joins the parts in order), so the order of a composition series is set on
    the command line. A categorical x axis uses this order; do not sort it.
    """
    return list(dict.fromkeys(df["file"]))


# ── Colours ──────────────────────────────────────────────────────────────────

def color_map(names):
    """Map names onto the style's colour cycle, wrapping around.

    Callers pass ``files_in(df)``, so one trajectory has one colour in every figure.
    """
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return {n: cycle[i % len(cycle)] for i, n in enumerate(names)}


# ── Layout ───────────────────────────────────────────────────────────────────

def panel_grid(n, cfg, wspace=None):
    """A grid of ``n`` panels laid out by ``cfg`` (``ncols``, ``panel_w``, ``panel_h``,
    ``wspace``, ``hspace``), the cells past ``n`` hidden; returns ``(fig, flat_axes)``.

    ``ncols`` None puts every panel in one row. ``wspace`` overrides ``cfg["wspace"]``.
    """
    ncols = cfg["ncols"] or n
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(cfg["panel_w"] * ncols, cfg["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=cfg["wspace"] if wspace is None else wspace,
                        hspace=cfg["hspace"])
    flat = axes.ravel()
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig, flat


# ── Display names and export ─────────────────────────────────────────────────

LABELS_FILE = "labels.csv"


class Context:
    """The shared state of one figure: display names, the compositions shown and their order,
    colours, the rows to export.

    ``files`` are the compositions in all inputs of this run (order of appearance), ``groups``
    the group names of several csv side by side. Both go into ``labels.csv``; only ``files``
    take part in the colour assignment.
    """

    def __init__(self, files, outdir=None, groups=()):
        self.path = (Path(outdir) if outdir else Path.cwd()) / LABELS_FILE
        table = load_labels(self.path, [*files, *groups])
        self._label = {k: lab for k, lab, _ in table}
        shown = {k for k, _, show in table if show}
        present = set(files)
        self.names = [k for k, _, _ in table if k in present and k in shown]
        if not self.names:
            raise ferro.FerroError(f"{self.path} sets show=0 for every composition of this run: "
                                   "nothing to draw")
        self._groups_shown = {g for g in groups if g in shown}
        self.colors = color_map(self.names)
        self.rows = []

    def shown(self, df):
        """The compositions in df with ``show=1``, in the row order of ``labels.csv``."""
        present = set(df["file"])
        return [n for n in self.names if n in present]

    def group_shown(self, group):
        return group in self._groups_shown

    def label(self, key):
        return self._label.get(key, key)

    def record(self, **row):
        """Register one row of export data. A ``label`` column follows the ``file`` column."""
        out = {}
        for k, v in row.items():
            out[k] = v
            if k == "file":
                out["label"] = self.label(v)
        self.rows.append(out)

    def record_curve(self, panel, file, quantity, x, y):
        """Register a curve: one row per point, columns ``panel, file, label, quantity, x, y``."""
        lab = self.label(file)
        self.rows += [{"panel": panel, "file": file, "label": lab, "quantity": quantity,
                       "x": float(a), "y": float(b)} for a, b in zip(x, y)]


def context(frames, outdir=None, groups=()):
    """frames: [(path, df)]. The compositions are every ``file`` in all products, in order of
    appearance."""
    files = list(dict.fromkeys(f for _, df in frames for f in files_in(df)))
    return Context(files, outdir, groups)


def load_labels(path, keys):
    """Read ``labels.csv``, append the keys new to this run, return [(key, label, show)].

    The file is written only when there are new keys; existing rows go back unchanged (order,
    edited labels, show).
    """
    table = []
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            for rec in csv.DictReader(fh):
                show = str(rec.get("show", "1")).strip()
                if show not in ("0", "1"):
                    raise ferro.FerroError(f"{path}: key '{rec['key']}' has show '{show}'; "
                                           "it must be 0 or 1")
                table.append((rec["key"], rec["label"], show == "1"))
    known = {k for k, _, _ in table}
    new = [k for k in dict.fromkeys(keys) if k not in known]
    if new:
        first = not path.exists()
        table += [(k, k, True) for k in new]
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["key", "label", "show"])
            w.writerows((k, lab, int(show)) for k, lab, show in table)
        if first:
            note(path, "written; edit label / show / row order and run again to apply them")
        else:
            note(path, f"{len(new)} new key(s) appended: {', '.join(new)}")
    return table


def save(fig, name, outdir=None, ctx=None):
    """Write the pdf (the product) and the png (for a look), and ``<name>_data.csv`` when rows
    were registered. Returns the paths written.

    ``bbox_inches="tight"`` comes with the ``science`` style and is not repeated here: an
    explicit figsize would otherwise lose its margins and be distorted.
    """
    out = Path(outdir) if outdir else Path.cwd()
    out.mkdir(parents=True, exist_ok=True)
    pdf, png = out / f"{name}.pdf", out / f"{name}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=PNG_DPI)
    written = [pdf, png]
    if ctx is not None and ctx.rows:
        data = out / f"{name}_data.csv"
        pd.DataFrame(ctx.rows).to_csv(data, index=False)
        written.append(data)
    plt.close(fig)
    return written
