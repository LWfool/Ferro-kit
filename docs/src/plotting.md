# Plotting

ferro does not draw figures. Every analysis writes a long-table csv (one row per
point, a `file` column naming the input), and figures are made from those files in
Python.

## Publication scripts

`scripts/` ships four matplotlib scripts built on a shared style layer
(`ferroplot.py`: SciencePlots `science` + `vibrant`, LaTeX text, PDF output). They need
`matplotlib`, `pandas`, `scienceplots` and a working LaTeX installation.

| Script | Reads | Figure |
|---|---|---|
| `plot_gr.py` | `gr_<pair>.csv` | one panel per csv; g(r) solid on the left axis, CN(r) dashed on the right |
| `plot_sq.py` | `sq.csv` | the X-ray and neutron totals, one row per csv |
| `plot_angle.py` | `angle_<triplet>.csv` | one panel per csv, one curve per input |
| `plot_net.py` | `network_*.csv` | 100 % stacked bars over composition |

```bash
python scripts/plot_gr.py gr_P-O.csv gr_Al-O.csv --outdir figs
python scripts/plot_sq.py sq.csv --outdir figs
python scripts/plot_angle.py angle_O-P-O.csv --outdir figs
python scripts/plot_net.py qn network_qn.csv --outdir figs
```

Colours follow the `file` column, so one input keeps one colour across every panel
and figure.

There is no script for MSD yet. The long table already works with one line of seaborn:

```python
sns.lineplot(data=pd.read_csv("msd_all.csv", comment="#"), x="time", y="msd", hue="file")
```

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
