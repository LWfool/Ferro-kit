# Merging Datasets

`ferro dataset merge` combines DeePMD systems of the same composition into one,
optionally shuffling them, and cuts the result into sets.

```bash
ferro dataset merge -i data/*.train -o merged
ferro dataset merge -i data/*.train -o merged --seed 42
ferro dataset merge -i data/*.train -o merged --mode by-source
ferro dataset merge -i clean/* -o merged --set-size 250 --suffix .train
```

## Grouping is by content, not by name

`init.011` tells you nothing reliable about what is inside it. Systems are
grouped by their **per-atom element sequence**, so two of them merge exactly
when they hold the same atoms — and systems of different composition separate
automatically, however they are named.

The output directory is `<natoms>_<formula>`:

```
Al2ZnO4_1.000_0300K.train  ┐
init.011                   ├─→  7_Al2O4Zn.train
whatever                   ┘
```

Subscripts are **actual counts, not reduced** — `Al96O192Zn48` is 336 atoms —
and the atom-count prefix makes `ls` group systems of equal size together.

The suffix is inherited when every input in the group shares one of `.train`,
`.test`, `.valid`; `--suffix` overrides, and inputs without a common suffix
produce a bare name.

## Atom order is canonicalised

Systems of one composition may still list their atoms in different orders and
carry different `type_map` orders. Merging sorts every system into one canonical
order — **(Z, symbol)** — and permutes the per-atom arrays along with them:

| Data | On reorder |
|---|---|
| `coord`, `force` | permuted with the atoms |
| `box`, `energy`, `virial` | untouched — independent of atom numbering |

A DP model is invariant under atom renumbering, so this changes notation rather
than physics; it is the same thing dpdata does with `sort_atom_names` plus
`sort_atom_types`. `type.raw` and `type_map.raw` are rewritten to match.

**ferro sorts by (Z, symbol) where dpdata sorts alphabetically.** Both are
self-describing through `type_map.raw`, and (Z, symbol) is what `collect` writes
and what `ferro traj gr` groups by — internal consistency wins over matching
dpdata's choice.

## The two modes

### `--mode shuffle` (default)

Everything of one composition is concatenated, shuffled with `--seed` (default
666, so a run without the flag is still reproducible), then cut into sets. Every
set holds a mix of whatever went in — temperatures, compressions, sources.

### `--mode by-source`

No mixing, no shuffling. Every source is cut into sets on its own and **the set
boundaries fall exactly on source edges**, so each `set.NNN` holds frames from a
single condition. The mapping is written to `sets_source.txt`:

```
# set  source
set.000  data/run300K.train
...
set.005  data/run500K.train
```

Merging 2000 + 500 frames at `--set-size 400` gives `400×5 + 250×2` — five sets
from the first source, two from the second, and no set straddling the join.

Pick `by-source` when the sources are already shuffled internally and you want a
validation split that is a clean hold-out of one condition; pick `shuffle` when
you want every set to be statistically like every other.

In both modes a remainder is spread across the sets rather than left as a stub,
because a set of a dozen frames is useless as a validation split.

## The pipeline

```
ferro dataset collect  -i md.out    -o raw      # AIMD  -> dataset
ferro dataset filter   -i raw       -o clean    # drop bad frames
ferro dataset merge    -i clean/*   -o merged   # combine + shuffle + resize
```

Each step reads and writes the same DeePMD directory format, and none of them
modifies its input.
