# Plotting

ferro does not draw figures. Every analysis writes a long-table csv (one row per
point, a `file` column naming the input), and figures are made from those files in
Python.

## Publication scripts

`scripts/` ships five matplotlib scripts built on a shared style layer
(`ferroplot.py`: SciencePlots `science` + `vibrant`, LaTeX text, PDF output). They need
`matplotlib`, `pandas`, `scienceplots` and a working LaTeX installation.

| Script | Reads | Figure |
|---|---|---|
| `plot_gr.py` | `gr_<pair>.csv` | one panel per csv; g(r) solid on the left axis, CN(r) dashed on the right |
| `plot_sq.py` | `sq.csv` | the X-ray and neutron totals, one row per csv |
| `plot_angle.py` | `angle_<triplet>.csv` | one panel per csv, one curve per input |
| `plot_msd.py` | `msd_<elements>.csv` | one panel per csv; fitted line over the fit window, $D \pm d_{err}$ in the legend; `--components` adds x/y/z, `--loglog` adds a slope-1 guide |
| `plot_net.py` | `network_*.csv` | four kinds, see [below](#network-figures) |

```bash
python scripts/plot_gr.py gr_P-O.csv gr_Al-O.csv --outdir figs
python scripts/plot_sq.py sq.csv --outdir figs
python scripts/plot_angle.py angle_O-P-O.csv --outdir figs
python scripts/plot_msd.py msd_all.csv --loglog --outdir figs
python scripts/plot_net.py qn network_qn.csv --outdir figs
python scripts/plot_net.py linkmap network_linkage.csv --outdir figs
```

Colours follow the `file` column, so one input keeps one colour across every panel
and figure.

### Display names: `labels.csv`

The first run writes `labels.csv` into the output directory (`--outdir`, or the
current directory). **Every script reads the same file**, so a composition is
renamed once for the whole paper.

```csv
key,label,show
43Z43P15A_NPT_5,$x=0.15$,1
70Z30P00A_NVT_5,$x=0$,1
```

- `label` is what the figure shows (LaTeX allowed); it starts equal to `key`
- `show=0` hides that input
- **row order is display order**: x axis, legend and colour assignment
- `key` is an input name (`file` column) or, when several csv are plotted side by
  side, a group name (the `-s` suffix: `CMD`, `MLMD`)

Later runs never touch existing rows. Keys seen for the first time are appended at
the end and reported with a `Note:` line — a misspelt key shows up there. Rows for
inputs that are absent from the current run are ignored.

### Exported data: `<stem>_data.csv`

Every figure is written together with a long table of exactly the numbers drawn
(after smoothing and filtering), for redrawing in a journal's style or next to
experimental data. It carries both the raw `file` and the mapped `label`.

| Figure type | Columns |
|---|---|
| curves (gr, sq, angle, msd) | `panel, file, label, quantity, x, y` |
| stacked bars (`plot_net.py qn / cn / bridge`) | `panel, group, file, label, category, percent, count` |
| heat maps (`plot_net.py linkmap`) | `panel, group, file, label, row, col, percent, count` |

```python
df = pd.read_csv("figs/net_bridge_data.csv")
df.pivot_table(index="label", columns="category", values="percent")
```

## Network figures

`plot_net.py <kind>` reads the products of `ferro net`.

| Kind | Reads | Figure |
|---|---|---|
| `qn` | `network_qn.csv` | one panel per Qn former; 100 % stacked bars over composition, bands $Q^0 \ldots Q^4$ |
| `cn` | `network_coordination.csv` | one panel per element with more than one coordination number |
| `bridge` | `network_linkage.csv` | left: the four ligand species `O_f` / `O_n` / `O_b` / `O_t`; right: the share of each former–former linkage (Al–O–Al, Al–O–P, P–O–P …) |
| `linkmap` | `network_linkage.csv` | heat maps; rows are compositions, columns are former pairs (3 for two formers, 6 for three), cells are the states at both ends |

`bridge` finds the ligand-species numbers in the `network_composition` csv beside the
linkage file, so only the linkage file is passed.

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
  to show; both are skipped.

`plot_msd.py` does not fit anything: the slope, intercept, $D$ and its error come from the
`[inputs]` list ferro writes, so there is one implementation of the fit.

## Why ferro does not plot

Up to `0.3.2` the trajectory commands took `--plot` and wrote a PNG through the Rust
`plotters` crate. It was removed in `0.3.3`:

- **It broke builds on Linux.** Text rendering needs the system `fontconfig` and
  `freetype` development libraries, which compute clusters often lack. On macOS the
  system font stack hid the problem.
- **It could not reach publication quality.** `plotters` has no PDF backend, and
  vector output through SVG → PDF would add about 40 more crates. The PNG was frozen
  at self-check level while the real figures were already made in Python.
- **It cost more than it gave.** It pulled in over 30 of the roughly 130 crates in
  the build, for four figure types that matplotlib does better.
