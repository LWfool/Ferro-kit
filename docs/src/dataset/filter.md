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
all frames  ->  |F|max  ->  |σ|max  ->  [start:end:stride|number]
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

Three tables are printed (not written — the dataset is the product; these
statistics exist to help you pick thresholds and change every run):

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

## Output

Survivors are written in their original order (**not shuffled**) and split into
sets of `--set-size` frames (default 400, `0` keeps one set). The remainder is
spread across the sets rather than left as a stub: 410 frames at `--set-size
400` give 205 + 205, not 400 + 10, because a set of ten frames is useless as a
validation split.

Writing into an existing non-empty directory requires `--overwrite`.

## Not implemented yet

- **Geometric criteria** — minimum O–O distance and Al⁶ coordination selection
  (`--oo-min`, `-r` in the reference Python implementation). They need periodic
  neighbour searching and their own diagnostic tables.
- **Shuffling and merging** — `ferro dataset merge`.
- **Extra keys.** `atom_ener.npy` and friends have no home in `Frame`; they are
  **reported** on read and would be lost on write. The reader names them rather
  than dropping them quietly.
