# Filtering a Dataset

`ferro dataset filter` drops low-quality frames from DeePMD system directories
and writes the survivors as a new dataset. **The input is never modified.**

```bash
ferro dataset filter -i raw                      # look, write nothing
ferro dataset filter -i raw -o clean
ferro dataset filter -i raw -o clean -f 15 -s 8
ferro dataset filter -i raw -o clean -N 500 --set-size 250
```

`-i` takes system directories or a directory holding them — a directory
containing `type.raw` is taken as a system and not descended into. `-o` is the
output **root**; each system is rebuilt under its path relative to `-i`, so
`raw/500K/a.train` becomes `clean/500K/a.train` and two systems with the same
name under different parents cannot overwrite each other.

**Omitting `-o` is the read-only mode.** It reports and writes nothing, which is
how you choose the thresholds before committing to them.

## The funnel

```
all frames -> |F|max -> |σ|max -> min d(O-O) -> Al6 -> [start:end:stride|number] -> shuffle
```

| Criterion | Flag | Default | Quantity |
|---|---|---|---|
| force | `-f` | 20 eV/Å | largest force **vector magnitude** over the atoms of a frame |
| stress | `-s` | 10 GPa | largest **absolute value among the 9 stress components** |

A threshold of `0` switches that criterion off. An explicit zero says "do not
judge", which no small positive number can express.

The stress criterion catches both diagonal and shear outliers, and in practice
it flags different frames than the force one — they are two independent checks,
not one plus a redundancy. Whether that holds for *your* data is exactly what
the report answers.

Stress is recovered as `virial / volume` with the volume from `|det(box)|`; a
DeePMD system carries no `volume.npy`. `-s` is given in GPa for convenience and
converted internally to eV/Å³.

## The range counts survivors

`--start` / `--end` / `--stride` / `-N` operate on the frames that **passed the
quality criteria**, not on original frame numbers. `--start 10` means "the 10th
surviving frame" — the only reading that stays meaningful after an unknown
number of frames were dropped.

`--end` is 0-based and **inclusive**, matching `ferro convert`. `--stride` is a
spacing, `-N` is a total; they are mutually exclusive, and `-N` always takes
both endpoints.

The report names original indices throughout, so kept frames remain traceable
to the input dataset.

## Reading the report

Three tables account for the selection itself:

| Table | What it says |
|---|---|
| `[funnel]` | how many frames each step left, in execution order |
| `[criteria]` | per criterion: frames flagged, and how many of those **no other criterion flagged** |
| `[overlap]` | frames flagged by both members of each criterion pair |

**The exclusive count is the one that matters.** The funnel alone cannot tell a
useful criterion from a redundant one, because each step only reports on what
the previous step left — a criterion that merely re-catches another one's frames
still looks productive there. An exclusive count near zero means the criterion
is not earning its place.

A real example, on a 2000-frame AIMD dataset with `-f 3 -s 0.2`:

```
[criteria]
  criterion  flagged  exclusive
  force      1978     0
  stress     2000     22
[overlap]
  a      b       both
  force  stress  1978
```

Every frame `force` flagged was also flagged by `stress`. At these thresholds
the force criterion contributes nothing, and the funnel — which would show
`2000 → 22 → 0` — hides that completely.

## Frame order

Kept frames stay in **trajectory order** unless `--shuffle` is given, and the
shuffle runs **last** — after every criterion and after `--stride` / `--number`.

That order is not arbitrary. "Take every 3rd frame" says nothing about a
shuffled sequence, so sampling has to see time order; and once shuffled, time
order cannot be recovered. Sorting the kept frames of a `-N 500 --shuffle` run
back into order gives exactly the same `[0, 4, 8, … 1999]` as the run without
it — the shuffle changes the write order, not which frames were chosen.

`--seed` (default 666) makes it reproducible, and is rejected without
`--shuffle` rather than silently ignored.

Shuffling here and shuffling in `ferro dataset merge` are **alternatives, not a
sequence**: shuffle in `merge` when sets should mix several sources, shuffle
here when this dataset goes to a trainer as it is.

## Output

Survivors are split into sets of `--set-size` frames (default 400, `0` keeps one set). The remainder is
spread across the sets rather than left as a stub: 410 frames at `--set-size
400` give 205 + 205, not 400 + 10, because a set of ten frames is useless as a
validation split.

Writing into an existing non-empty directory requires `--overwrite`.

### The report files

With `-o`, all seven tables are also written as CSV, flat in the output root:

```
clean/
├── 49Z49P02A_0.950_0300K/     the filtered systems, one per input
├── 49Z49P02A_0.950_1000K/
├── filter_funnel.csv          ┐
├── filter_criteria.csv        │ the selection
├── filter_overlap.csv         ┘
├── filter_min_oo.csv          ┐
├── filter_al6.csv             │ the diagnostics
├── filter_al_cn.csv           │
└── filter_rcut_scan.csv       ┘
```

They go through the same writer as every other ferro product, so each carries a
`#` header with the shared parameters and an `[inputs]` list of the systems, and
`pandas.read_csv(comment="#")` reads them directly. Several systems stack into
one file with a `system` column holding the path **relative to `-i`** — nested
systems `a/md` and `b/md` share a leaf name and would otherwise be
indistinguishable once stacked.

They sit flat in the root rather than in a `report/` subdirectory on purpose: a
later `ferro dataset merge -i clean/*` keeps only directories, so the CSVs are
filtered out by themselves, while a `report/` directory would be picked up as a
system candidate.

**Without `-o` nothing is written at all** — every table is printed instead, the
four diagnostics included. With `-o` the diagnostics are written but not printed,
since a few dozen lines of them would scroll the destination path away.

## The two geometric criteria

Force and stress are the first-line rules, on by default. These two are optional
and independent — both, either, or neither.

### `--al6` — keep frames containing a 6-coordinated Al

"Keep frames containing at least one Al⁶" and "drop frames containing none" are
the same rule; ferro expresses it as the latter so all four criteria share one
sense and the cross-tabulation stays a single table.

Coordination comes from the **same classifier `ferro net` uses**
(`classify_frame`), so the number means the same thing in both commands. Only
the Al–O cutoff is needed: Al is not in the default Qn element set `{B,P,Si}`,
so it takes the "non-Qn former" branch where the label's digit *is* the
coordination number.

The cutoff can be derived. Bare `--al6` takes **the first minimum of the Al–O
g(r) past its first peak** — the outer edge of the first coordination shell —
computed per system, because compositions differ and so do shell positions. On
the reference system that gives 2.45 Å against the 2.4 Å normally used by hand.
The value is always printed, and the mean over systems reported at the end: a
cutoff that decides which frames die must not be an invisible number.

### `--oo-min` — drop frames with a too-short O–O contact

**This one has no automatic form and always needs an explicit value** (bare
`--oo-min` uses 2.0 Å). The reason is structural, not an omission:

| pair | first peak | first minimum past it | depth there |
|---|---|---|---|
| Al–O | 1.81 Å | **2.41 Å** | 0.23 |
| P–O | 1.51 Å | 2.25 Å | 0.08 |
| O–O | 2.55 Å | 3.77 Å | **0.72** |

O–O does not bond, so its g(r) has no coordination shell — the "first minimum"
is a shallow feature almost 4 Å out. And healthy frames have their smallest O–O
distance *below* the first O–O peak (1.77–2.14 Å against a peak at 2.55 Å). The
threshold you actually want comes from the RDF of **broken** data: the trough
between the collapse peak and the normal one. That trough does not exist in data
that is still good, which is why no amount of analysis of a healthy dataset can
produce it.

## The diagnostics

Four more tables help pick those two values. They are always computed: on
1110 frames of 302 atoms they cost no measurable wall time (0.23 s with them
against 0.28 s without — the work is spread over every core), and a table that
exists only in read-only mode is a table that never reaches an archive.

**The cutoff scan is the important one.** On one reference system it goes

```
rcut  frames_pct  per_frame
2.15  0.9         0.01
2.45  13.1        0.14
2.75  41.4        0.46
```

and on another it is a flat 100% from 2.1 to 2.6 Å. A steep column means the
selection is decided by the cutoff rather than by the structure. **The same
table says opposite things about the two systems** — which is exactly why it is
worth printing rather than baking a number into the code.

The `min d(O-O)` table reports the **distribution** — quantiles and a histogram
— rather than a count below the threshold. A count cannot tell an outlier tail
from a smooth spread, and only the first is worth filtering away. Note also that
this quantity is an *extreme-value* statistic: its location shifts with the
number of O atoms in the cell (64 O gave a mean of 2.40 Å, 201 O gave 2.00 Å for
the same kind of glass), so a threshold carried over from another system means
very little.

## Not implemented yet

- **Multi-layer periodic images.** Cutoffs are checked against the
  minimum-image bound and rejected beyond it, rather than scanning further image
  shells. Not a limitation for cells much larger than the cutoff.
- **Extra keys.** `atom_ener.npy` and friends have no home in `Frame`; they are
  **reported** on read and would be lost on write. The reader names them rather
  than dropping them quietly.
