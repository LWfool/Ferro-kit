# Plotting

ferro does not draw figures. Every analysis writes a long-table csv (one row per
point, a `file` column naming the input), and figures are made from those files in
Python.

Publication figures of `ferro traj gr / sq / angle / msd` and `ferro net` are drawn by
the Python side, `ferro-py plot` (SciencePlots style, LaTeX text, PDF output, with the
numbers drawn exported next to each figure):

```bash
ferro-py plot gr -i gr_P-O.csv gr_Al-O.csv -o figs
ferro-py plot net qn -i network_qn.csv -o figs
```

Its manual is `ferro-py doc plot`.

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
