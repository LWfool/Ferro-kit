# Building Structures

`ferro build` makes a new structure from existing ones, or from nothing (`box`).

```bash
ferro build supercell -i quartz.cif --dims 2 2 3
ferro build supercell -i 'cifs/*.cif' --min-atoms 200 --min-length 12 -o big --mkdir
ferro build vacuum -i slab.vasp --thickness 15
ferro build merge -i substrate.vasp water.vasp --gap 2.5
ferro build box -c ZnO:43 -c P2O5:43 -c Al2O3:15 --density 3.2
```

## Shared conventions

| Flag | Default | Description |
|---|---|---|
| `-o <DIR>` | current dir | output directory; a missing one is asked about, or created with `--mkdir` |
| `-s <SUFFIX>` | | batch tag, the last segment of every product name |
| `--mkdir` | | create `-o` without asking (required when there is no terminal) |
| `--format <EXT>` | the input's | output format as an extension: `vasp`, `cif`, `extxyz`, `data`, … |

Products are named like every other ferro product, `<op>_<input>[_<suffix>].<ext>`:
`supercell_quartz.cif`, `vacuum_slab_s1.vasp`. `merge` names both inputs
(`merge_<A>_<B>`), and `box`, which has none, is `box[_<suffix>].vasp`.

- **The format follows the input.** A `POSCAR` / `CONTCAR` input gives `.vasp`. An
  input ferro only reads (OUTCAR, a program output, a QE input) needs `--format`,
  and the run stops before anything is read if it is missing.
- **A single-structure format gets one file per frame**, as with `convert`:
  a two-frame trajectory written as POSCAR is `supercell_md_0000.vasp`,
  `supercell_md_0001.vasp`.
- **Atoms are written in (Z, symbol) order**, one block per element, as `convert` writes them.
- **Energy, forces, stress, velocities, temperature and step are dropped.** They
  belong to a calculation, and the new structure has not been calculated.
  Labels, masses, magnetic moments and charges are kept.
- `supercell` and `vacuum` take several inputs (`-i` takes globs; quote them). An input
  that fails is skipped and the exit code is 1, as in every batch command.

## `ferro build supercell`

Replicates a cell $n_x \times n_y \times n_z$ times along $a$, $b$, $c$.

| Flag | Default | Description |
|---|---|---|
| `-i <FILE>...` | | input structure(s) |
| `--dims <NX> <NY> <NZ>` | | the replicas; excludes the two below |
| `--min-atoms <N>` | | the supercell holds at least $N$ atoms |
| `--min-length <LEN>` | | every interplanar spacing of the supercell is at least LEN Å |
| `--dry-run` | | print the replicas and atom counts, write nothing |

Without `--dims`, the replicas are the smallest that meet both conditions:
$n_i = \lceil L / d_i \rceil$ for each axis, with $d_i$ the interplanar spacing of the
input cell, then the axis with the thinnest supercell grows by one until the atom
count reaches `--min-atoms`. Growing the thinnest axis keeps the box as close to
isotropic as it can be.

**`--min-length` is a spacing, not a vector length.** The spacing $d_i$ is the
distance between the lattice planes spanned by the other two vectors — the
thickness of the cell along that axis. It is what bounds a minimum-image cutoff,
$r_\text{cut} \le d_\text{min} / 2$, so "every spacing at least 12 Å" is the condition
for a 6 Å cutoff. On an orthogonal cell the two are the same; on a tilted cell the
spacing is shorter, and sizing by vector length gives a cell that is too thin.

The replicas are worked out from the **first frame** and applied to every frame, so
all frames have the same number of atoms. When the input is already large enough,
`1x1x1` is printed and the copy is written all the same: a script that expects the
product finds it.

The supercelled trajectory is held in memory at once, so it is capped at $10^8$ atoms
in total (frames times atoms per frame); `--min-atoms` above that fails before any
input is read. A long trajectory over the cap is refused before anything is built:
supercell a single frame instead.

The line printed for each input shows the replicas, how they were chosen (`given`
or `auto`), the atom counts and the supercell's three spacings:

```
[1/1] POSCAR
  4x3x3 (auto): 3 -> 108 atoms, spacings 38.509 35.909 42.000 Å -> supercell_POSCAR.vasp
```

## `ferro build vacuum`

Adds a vacuum layer along one lattice vector.

| Flag | Default | Description |
|---|---|---|
| `-i <FILE>...` | | input structure(s) |
| `--axis <X>` | z | `x` = $a$, `y` = $b$, `z` = $c$ |
| `--thickness <LEN>` | | vacuum thickness in Å, greater than 0 |

The chosen vector is lengthened, keeping its direction, until the interplanar
spacing along it has grown by `--thickness`: the gap between periodic images is the
thickness asked for, even when the vector is tilted. The other two vectors do not
change.

**Atoms are not moved.** A slab that crosses the cell boundary along `--axis` is cut
in two by the gap; move it into the middle of the cell first.

## `ferro build merge`

Stacks two structures along one lattice vector, with a gap between them.

| Flag | Default | Description |
|---|---|---|
| `-i <A> <B>` | | the two structures, bottom first, one frame each |
| `--axis <X>` | z | the lattice vector to stack along |
| `--gap <LEN>` | | gap between the two blocks in Å, measured along the interface normal; 0 is allowed |

**Neither block is strained.** Both are moved rigidly, and in the interface plane
each is centred on the merged cell. The merged cell keeps the directions of A's
lattice vectors; each in-plane vector is as long as the longer of the two blocks',
and the stacking vector is stretched until the height along the interface normal is
$h_A + \text{gap} + h_B$, with $h$ each block's interplanar spacing. This makes a padded
starting box, not a coherent interface: ASE `stack` and pymatgen
`Interface.from_slabs` strain B onto A's lattice instead, which is impossible when
the in-plane sizes differ a lot.

- The two interfaces must be parallel (normals within 1°): a gap along the normal
  means nothing otherwise.
- Each block is first folded into its own cell along its periodic axes.
- A trajectory is refused rather than reduced to its first frame. Take the frame you
  want with `ferro convert -i md.dump -o DIR --format extxyz --start K --end K`.
- The output format is A's unless `--format` says otherwise.

## `ferro build box`

Packs molecules at random into a cubic box of a given density.

| Flag | Default | Description |
|---|---|---|
| `-c <NAME:COUNT>` | | one component; repeat for each |
| `--density <RHO>` | | mass density in g/cm³ |
| `--min-dist <LEN>` | 1.5 | closest allowed distance between two atoms, in Å |
| `--relax-steps <N>` | 500 | relaxation steps |
| `--seed <N>` | 666 | random seed |
| `--format <EXT>` | vasp | output format |

A component is a formula — nested brackets allowed, `Ca3(PO4)2`, `K4[Fe(CN)6]` — or
one of the known solvent names (`water`, `methanol`, `ethanol`, `acetone`, `dmso`,
`dmf`, `thf`, `benzene`, `toluene`, …), whose tabulated molar mass is used. Quote a
formula with brackets for the shell. Hydrate notation ($CuSO_4 \cdot 5H_2O$) is not
accepted: list the water as its own component.

The box length is $L = (m / \rho)^{1/3}$ from the total mass $m$. The atoms are placed
at random, then a soft-core relaxation pushes apart every pair closer than
`--min-dist`. If such pairs remain after `--relax-steps`, the run fails rather than
write a box with overlapping atoms: raise `--relax-steps`, or lower `--density`.

**Molecules are not kept whole.** Each atom is placed on its own, so this is a
starting structure for a melt-quench run, not a molecular liquid. The same seed
gives the same box, atom for atom.
