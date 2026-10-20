# ferro-py plot: Publication Figures

```bash
ferro-py plot gr    -i gr_P-O.csv gr_Al-O.csv -o figs
ferro-py plot sq    -i sq.csv -o figs
ferro-py plot angle -i angle_O-P-O.csv -o figs
ferro-py plot msd   -i msd_all.csv --loglog -o figs -s loglog
ferro-py plot net qn      -i network_qn.csv -o figs
ferro-py plot net linkmap -i network_linkage.csv -o figs
```

ferro writes csv only; `ferro-py plot` turns those products into figures for a paper:
SciencePlots `science` + `vibrant` style, LaTeX text, a PDF with editable text. It needs
the optional dependencies and a working LaTeX installation:

```bash
pip install 'ferro[plot]'      # matplotlib, pandas, scienceplots
```

| Figure | Reads | Draws |
|---|---|---|
| `gr` | `gr_<pair>.csv` | one panel per csv; $g(r)$ solid on the left axis, CN(r) dashed on the right |
| `sq` | `sq.csv` | the neutron and X-ray totals, one row of two panels per csv |
| `angle` | `angle_<triplet>.csv` | one panel per csv, $P(\theta)$ |
| `msd` | `msd_<elements>.csv` | one panel per csv; fitted line over the fit window, $D \pm d_{err}$ in the legend; `--components` adds x/y/z, `--loglog` adds a slope-1 guide |
| `net` | `network_*.csv` | four tables, see [below](#network-figures) |

Every figure takes the same options: `-i` the csv files (quoted globs allowed), `-o` the
output directory (the current one by default; a missing one is asked about, `--mkdir`
when not interactive) and `-s` the last segment of the name. A figure named `<name>`
is written as `<name>.pdf`, `<name>.png` (a preview) and `<name>_data.csv`; `<name>` is
`gr`, `sq`, `angle`, `msd` or `net_<table>`, followed by `_<-s>` when given.

Colours follow the `file` column, so one input keeps one colour across every panel
and figure. A gr or angle csv holding several pairs or triplets (`gr_all.csv`) is refused:
one panel cannot say which one it shows.

`msd` does not fit anything: the slope, intercept, $D$ and its error come from the
`[inputs]` list ferro writes, so there is one implementation of the fit.

## Display names: `labels.csv`

The first run writes `labels.csv` into the output directory. **Every figure reads the
same file**, so a composition is renamed once for the whole paper.

```csv
key,label,show
43Z43P15A_NPT_5,$x=0.15$,1
70Z30P00A_NVT_5,$x=0$,1
```

- `label` is what the figure shows (LaTeX allowed); it starts equal to `key`
- `show=0` hides that input
- **row order is display order**: x axis, legend and colour assignment
- `key` is an input name (`file` column) or, when several csv are plotted side by
  side, a group name (the `-s` suffix of `ferro net`: `CMD`, `MLMD`)

Later runs never touch existing rows. Keys seen for the first time are appended at
the end and reported as a `FerroNote` — a misspelt key shows up there. Rows for inputs
that are absent from the current run are ignored.

## Exported data: `<name>_data.csv`

Every figure is written together with a long table of exactly the numbers drawn
(after smoothing and filtering), for redrawing in a journal's style or next to
experimental data. It carries both the raw `file` and the mapped `label`.

| Figure | Columns |
|---|---|
| curves (gr, sq, angle, msd) | `panel, file, label, quantity, x, y` |
| stacked bars (`net qn / cn / bridge`) | `panel, group, file, label, category, percent, count` |
| heat maps (`net linkmap`) | `panel, group, file, label, row, col, percent, count` |

```python
df = pd.read_csv("figs/net_bridge_data.csv")
df.pivot_table(index="label", columns="category", values="percent")
```

## Network figures

`ferro-py plot net <table>` reads the products of `ferro net`.

| Table | Reads | Figure |
|---|---|---|
| `qn` | `network_qn.csv` | one panel per Qn former; 100 % stacked bars over composition, bands $Q^0 \ldots Q^4$ |
| `cn` | `network_coordination.csv` | one panel per element with more than one coordination number (`--element` to choose) |
| `bridge` | `network_linkage.csv` | left: the four ligand species `O_f` / `O_n` / `O_b` / `O_t`; right: the share of each former–former linkage (Al–O–Al, Al–O–P, P–O–P …) |
| `linkmap` | `network_linkage.csv` | heat maps; rows are compositions, columns are former pairs (3 for two formers, 6 for three), cells are the states at both ends |

`bridge` finds the ligand-species numbers in the `network_composition` csv beside the
linkage file, so only the linkage file is passed. Several csv of one table (runs tagged
`-s CMD`, `-s MLMD`) stand side by side, one bar per run on each composition.

What the numbers mean:

- **Shares come from summed counts**, not from the per-frame `fraction` column: after
  rows are merged the mean of per-frame fractions cannot be recovered, and `sd` does
  not add. No error bars are drawn.
- **The linkage share counts only bridges between formers**; non-bridging ligands are
  not in the denominator.
- **Triclusters count as connections.** `ferro net` expands a ligand shared by $k$
  formers into $C(k,2)$ pairs, and both figures use all of them, as $Q^n$ does. How
  many ligands are triclusters is shown by the `O_t` band.
- **Heat-map axes**: a Qn former is placed by $n$, its number of homonuclear
  connections; any other former is placed by its coordination number. In the Al–P
  map, $Q^0$ therefore means a P bonded to Al but with no P–O–P. Which elements are
  Qn formers is read from the `mean_qn` column of `[inputs]`.
- **Each map sums to 100 %** on its own; the colour scale is shared down a column so
  compositions compare directly. Maps of one element (P–O–P) show the upper
  triangle only, because a bridge has no direction and is stored once.
- **Grey cells**: a light-grey cell never occurs in any composition; a cell reading
  `0` occurs elsewhere but not here; an all-grey map marked *absent* means that
  composition has no such linkage.
- With a single former, the linkage panel is always 100 % and `linkmap` has nothing
  to show; both are skipped with a note.

## From Python

Each figure is a module of `ferro.plot` with a `plot` function returning the paths it
wrote; the style settings are the module's `CFG` dictionary.

```python
from ferro.plot import gr, net

gr.CFG["xlim"] = (0, 6)
gr.plot(["gr_P-O.csv", "gr_Al-O.csv"], outdir="figs")
net.plot(["network_coordination.csv"], "cn", outdir="figs", elements=["Al", "Zn"])
```

`import ferro` does not load matplotlib; importing a `ferro.plot` module without the
optional dependencies raises `ferro.FerroError`.
