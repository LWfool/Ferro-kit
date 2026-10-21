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

