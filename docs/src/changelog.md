# Changelog

What changed for *you* — the command lines you have to edit, the products you
have to regenerate. The reasoning behind each decision lives in the development notes,
not here.

Versions from `v0.3.0` onwards are listed. Anything produced by `0.2.x` or
earlier is incompatible with today's ferro in file names, column structure and
command line alike; regenerate it rather than migrating it.

---

## Unreleased — the `0.3.2` / `0.3.3` batch

Three rounds of breaking changes sit under the version number `0.3.3`, none of
them released yet. They will arrive together.

### `-o` is always a path

`--outdir` is gone. `-o` names a **directory** for every command whose run
writes several products (the 13 analysis commands, `map chg-sdf`, `bader`,
`dataset`), and the **output file** for `convert` and `job`, which write
exactly one. A batch tag goes to `-s` / `--suffix`.

| before | now |
|---|---|
| `ferro traj gr -i t.dump -a P -b O -o cmp` | `ferro traj gr -i t.dump -a P -b O -s cmp` |
| `ferro traj gr … --outdir scan` | `ferro traj gr … -o scan` |
| `ferro dataset filter -i raw --outdir clean` | `ferro dataset filter -i raw -o clean` |

`-o cmp` does not fail — it now means *write into the directory `cmp/`*. Check
every script that passed a tag to `-o`.

A missing directory is offered for creation on stderr; **a non-interactive run
(script, CI) must pass `--mkdir`** or it exits with an error.

### `ferro bader` writes beside its input

The three reports default to the **input file's own directory**, not the
current one, and are named after the input stem: `run1/CHGCAR` gives
`run1/CHGCAR_ACF.dat`. VASP calls every charge density `CHGCAR`, so the old
default made two runs overwrite each other. `-o` collects them elsewhere, `-s`
tags them.

### `dataset` product names

`collect` products carry `.db`, `filter` and `merge` products carry `.train`
unless `--ratio` splits them into `.train` / `.valid` / `.test`.

| before | now |
|---|---|
| `sets/run1/` | `sets/run1.db/` |
| `clean/md/` | `clean/md.train/` |

`merge` refuses a mix of split suffixes, and `.valid` / `.test` cannot be split
again. `collect`'s `--outdir` became `--output` (`-o` unchanged), and its `-o`
is optional: without it, products land **beside** each AIMD directory.

### One system per directory

`collect` groups by **directory**, not by file: the `.out` files of one
directory are the restart segments of one run and are reassembled into one
system. Two compositions in one directory is an error, not a dropped frame.

Atom order is **no longer** one of those differences. Every file is sorted into
canonical `(Z, symbol)` order as it is read, so two single points of one
material that merely list their atoms differently now collect into one system.
Products therefore have a different atom order than `0.3.1` produced — training
is unaffected (`type_map.raw` is self-describing, and DeePMD-kit re-sorts by
type on load anyway), but the files are not byte-identical.

### `collect --format`

```bash
ferro dataset collect -i 'sp/*.log' --format cp2k/sp
```

`cp2k/md` | `cp2k/sp` | `vasp/outcar` | `vasp/xml`. The banner is still read;
`--format` overrides it and a mismatch prints one `NOTE:` line. A CP2K file
with no `GLOBAL| Run type` line is now an error naming this flag, where it used
to be read as MD.

### LAMMPS triclinic dump boxes

The three box lines of a **triclinic** dump are written as the bounding box
(`xlo_bound`/`xhi_bound`), which is what the format specifies. Orthogonal cells
are byte-identical to before; triclinic ones were written too small, so OVITO,
ASE and LAMMPS `read_dump` all read the wrong cell. **Regenerate triclinic
dumps written by an earlier ferro.**

### extxyz stress sign

`stress=` follows the ASE convention (positive = tension). The nine numbers are
the negative of what `0.3.1` wrote. `virial=` is unchanged (eV, positive =
compression). A file carrying both is cross-checked and rejected if they
disagree; a 6-component Voigt vector is refused.

### New, nothing to migrate

`ferro doc <topic>` (the manual, compiled into the binary; rendered with tables and
formulas on a terminal, the plain source when redirected) · VASP `OUTCAR` and
`vasprun.xml` reading · CP2K single-point reading · `dataset collect --type
inspect` (diagnostics, no dataset) · `--type nep|extxyz` and `--ratio 8:1:1`
on `filter` and `merge` · `--qn` on `net`.

### Help pages are shorter

Every page is now five sections: what it does, the full parameter table, the
output layout, examples, and a pointer to the manual. Everything else moved
here into the manual. No command line changed.

---

## `v0.3.1` — machine-learning datasets (2026-08-26)

`ferro dataset collect | filter | merge`, DeePMD npy reading and writing, and
the CP2K MD reader. **Purely additive** — no existing behaviour changed.

---

## `v0.3.0` — glass network rework (2026-08-20)

Three rounds of breaking changes released together.

### `ferro net`: six tables, renamed and restructured

| before | now |
|---|---|
| `bridge` | `qn` |
| `partner` | `qn_partner` |
| `oxy` | `ligand_type` |
| `cn` | `coordination` |
| `mean` | removed |
| — | `composition` (new) |

`net qn` and `net type` are gone; `net` is a single leaf command. The number in
a non-Qn former's label is now its **coordination** number (`Al_4` = four
coordinated), not its bridge count.

### `ferro net`: `n` counts homopolar bridges only

The Qn convention follows the literature's Q^n_m: **`n` counts P–O–P bridges
only**, heteropolar ones go to the `m_<X>` columns of `qn_partner`, and the
total is `n + sum(m)`. **Every number in the `qn` table changed** — one
reference run moved from 40.4 % `P-Q3` to 0.27 %, and `mean_qn` from 2.40 to
0.95. Plots and comparisons built on the old numbers have to be redone.

### `ferro traj`: product names carry the selection

`gr.csv` became `gr_P-O.csv`, and `gr_all.csv` when nothing is selected. All
six commands changed, including the unselected case.

`traj sq` lost `-a` / `-b` / `-x` / `-y`: every pair is always written, because
the partials sum back to the weighted totals and keeping one pair would hide
exactly that. Filter columns in pandas instead.
