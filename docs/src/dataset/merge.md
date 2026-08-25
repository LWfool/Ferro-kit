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

No mixing, no shuffling. Each system is cut into sets **on its own**, so a set
never straddles two systems and every `set.NNN` stays traceable to one
condition. `--seed` does not apply. The mapping is written to
`sets_source.txt`:

```
# set  frames  source
set.000  400  data/run300K.train
...
set.004  400  data/run300K.train
set.005  250  data/run500K.train
set.006  250  data/run500K.train
```

Merging a 2000-frame and a 500-frame system at `--set-size 400` gives 5 + 2
sets. The point is that each set stays identifiable: if every input is one
condition — a temperature, a compression — then every set remains a clean
hold-out of that condition, and picking a validation split means picking a set.

Pick `shuffle` instead when you want every set to be statistically like every
other.

## The remainder is spread, in both modes

500 frames at `--set-size 400` give **250 + 250, not 400 + 100**. A lopsided
pair is worse both for training balance and for using a set as a validation
split — and a trailing set of a dozen frames is useless as either.

## The pipeline

```
ferro dataset collect  -i md.out    -o raw      # AIMD  -> dataset
ferro dataset filter   -i raw       -o clean    # drop bad frames
ferro dataset merge    -i clean/*   -o merged   # combine + shuffle + resize
```

Each step reads and writes the same DeePMD directory format, and none of them
modifies its input.
