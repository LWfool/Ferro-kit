# Collecting AIMD Output

`ferro dataset collect` turns *ab initio* MD output into
[DeePMD-kit](https://docs.deepmodeling.com/projects/deepmd/) system
directories — the starting point for training a machine-learning potential.

```bash
ferro dataset collect -i 'run*/*.out' -o data      # -> data/run1/, data/run2/
ferro dataset collect -i '*.out' -o sys            # -> sys/   (one directory in)
```

Three sources are read: **CP2K MD output**, **VASP `OUTCAR`** and **VASP
`vasprun.xml`**. Quantum ESPRESSO is still planned.

### Which reader gets the file

By **content**, not by name. The first lines of the file carry an unambiguous
banner — `CP2K|`, `vasp.6.4.2`, or an `<?xml` declaration — and naming cannot
be trusted to do this job: VASP writes `OUTCAR` with no extension at all,
people rename it to `run.outcar`, and `.out` is too generic to belong to any
one program. An unrecognised file is an error that names what *is* recognised.

> **One format per directory.** A real VASP run directory holds `OUTCAR` *and*
> `vasprun.xml`, and they record the same frames. Since `collect` treats the
> files of one directory as segments of one run, feeding it both would
> concatenate the same frames twice and silently double the dataset — the
> composition matches and both files parse, so nothing else would look wrong.
> Mixing formats within one directory is refused; narrow `-i` to one of them.

### CP2K vs VASP: what differs

| | CP2K | VASP OUTCAR | vasprun.xml |
|---|---|---|---|
| energy | `ENERGY\| Total FORCE_EVAL` | `free  energy   TOTEN` | last `e_fr_energy` of the calculation |
| convergence | `SCF run converged` | VASP's own `EDIFF is reached` | inferred: SCF steps < `NELM` |
| species | element column of the xyz block | `VRHFIN` × `ions per type` | `<atominfo>` |
| restarts | `MD_INI` blocks | ionic step counter going backwards | not detected |

The two VASP paths do **not** use the same convergence rule, so the same run
can drop a different number of frames depending on which file you point at.
The rule that applied is printed with the per-directory report rather than left
for you to guess.

`energy(sigma->0)` is deliberately *not* used: the forces VASP prints are the
derivatives of the free energy, so pairing them with the extrapolated energy
would give a model two halves of different functionals. dpdata makes the same
choice, which keeps datasets converted by either tool comparable.

## What CP2K must print

The run log must be self-contained — coordinates, forces and the stress tensor
all going to `__STD_OUT__` rather than to sibling files:

```
&MOTION
  &PRINT
    &TRAJECTORY
      &EACH
        MD 1
      &END EACH
      FILENAME __STD_OUT__
    &END TRAJECTORY
    &FORCES
      FILENAME __STD_OUT__
    &END FORCES
  &END PRINT
&END MOTION

&FORCE_EVAL
  STRESS_TENSOR ANALYTICAL
  &PRINT
    &STRESS_TENSOR
    &END STRESS_TENSOR
  &END PRINT
&END FORCE_EVAL
```

A restarted run whose logs were concatenated into one file is fine — restarts
are detected and reported, and CP2K does not reprint the initial configuration
on restart, so no duplicate frames arise.

Restart segments left as **separate files in one directory** are also fine; see
[One system per directory](#one-system-per-directory).

## Output layout

One system directory per input **directory**:

```
<outdir>/<name>/
├── type.raw          one 0-based integer per atom, indexing type_map.raw
├── type_map.raw      one element symbol per line, sorted by (Z, symbol)
└── set.000/
    ├── coord.npy     (nframes, natoms*3)   Å
    ├── box.npy       (nframes, 9)          Å, row-major lattice vectors
    ├── energy.npy    (nframes, 1)          eV
    ├── force.npy     (nframes, natoms*3)   eV/Å
    └── virial.npy    (nframes, 9)          eV
```

Every array on disk is **two-dimensional**: dpdata flattens with
`reshape([nframes, -1])` before saving, so the documented `nframes × natoms × 3`
is a logical shape, not the stored one. ferro follows the same convention, and
the files load unchanged with `numpy.load` or `dpdata.LabeledSystem`.

`<name>` is the path of the input directory **below the ancestor every input
shares**, kept nested rather than flattened with separators. A shared prefix
carries no distinguishing information by definition, so what is left after
stripping it is exactly what tells the systems apart:

| `-i` | products |
|---|---|
| `run*/*.out -o sets` | `sets/run1/`, `sets/run2/` |
| `/s/a/md/x.out /s/b/md/x.out -o sets` | `sets/a/md/`, `sets/b/md/` |
| `*.out -o sys` (one directory) | `sys/` itself |
| `a/total.out b/md/total.out -o sets` | `sets/a/`, `sets/b/md/` |

The file stem never enters the name — `total.out` and `PZA.out` in the same
directory produce the same system. With only one input directory the shared
ancestor is the whole path, `<name>` is empty, and the system is written into
`-o` itself: there is nothing to tell apart.

`-o` is **required**. The products are a directory tree, and a default of `.`
would scatter `.npy` files through whatever directory you happened to be in.
Writing into an existing non-empty directory needs `--overwrite`.

## One system per directory

The `.out` files sitting in one directory are the restart segments of one run,
so `collect` puts them back together into **one** system rather than one each.
That is the line between the two commands: `collect` reassembles the pieces of
**one** run, [`merge`](merge.md) combines **different** runs of the same
composition.

Files are ordered by their first `MD| Step number`, and each keeps its own
internal order. Sorting every frame globally would look more thorough, but a run
restarted without a checkpoint numbers its steps from zero again, and a global
sort would then interleave two real trajectories. The worst case here degrades
to "concatenate in file order", which is no worse than not sorting at all.

Overlapping frames are **not** removed. A restart re-runs at most the few steps
since the last checkpoint, and identical positions and velocities give identical
energies and forces, so the repeat neither biases nor dilutes the set. The step
span of every source file is printed so that premise stays checkable:

```
sets/run1  (2 file(s), 1021 frames)
  run1/a.out   steps 1-620     620 kept, 0 dropped
  run1/b.out   steps 500-900   401 kept, 0 dropped
```

Two files of **different composition** in one directory are an error naming both
files, not a frame-dropping event: a `type.raw` is written once per directory, so
the atom sequence must match throughout. Putting two systems in one directory is
a mistake of the person, not a problem with the data, and the two call for
completely different responses.

A file that fails to parse is skipped, the rest still become a system, and the
skipped files are listed again at the end with exit code 1. The second listing is
not redundant: the system directory looks perfectly normal while holding fewer
frames than you think.

### float64, not float32

dpdata defaults to `float32`; ferro writes `float64`. This directory is the head
of the pipeline — `filter` and `merge` read it back — and precision lost here
cannot be recovered downstream. Narrow to `float32` at the step that feeds the
training framework, not before. The cost is a factor of two in disk: for a
2000-frame, 112-atom run, 11 MB instead of 5.5 MB.

### One set, never split

`collect` always writes a single `set.000`. Splitting into `set.000`,
`set.001`… exists to give `merge --shuffle` its boundaries; nothing is shuffled
at collection time, and splitting early only makes `filter` read many small
files.

## Units and signs

| Quantity | In CP2K output | Stored |
|---|---|---|
| energy | `ENERGY\| Total FORCE_EVAL ( QS ) energy [hartree]` | eV |
| force | xyz block, **no unit printed** → atomic units | eV/Å |
| stress | `STRESS\| Analytical stress tensor [bar]` | eV/Å³ (→ virial in eV) |
| cell | the 12-field line after the force block | Å |

**Units are read from the text, not inferred from the version.** CP2K's
`STRESS_UNIT` is an *input* keyword, so one and the same binary can print bar,
GPa or atm. An unrecognised unit is an error — never a silent default. (dpdata
hard-codes GPa here; on a `bar` output that is wrong by five orders of
magnitude.) Forces are the one quantity CP2K prints with no unit at all, and are
taken as Hartree/Bohr.

**Sign**: the stress keeps the sign CP2K prints — *positive = compression* —
which is also the orientation DeePMD's virial uses, so `virial = stress × V`
with no flip. In ASE terms both equal `−V·σ_ASE`. VASP and Quantum ESPRESSO
print the same orientation as CP2K (ASE flips the sign when reading either).
GPUMD's `stress=` keyword, by contrast, uses the ASE convention and needs a
flip; its `virial=` keyword does not.

**The stress is the potential part only.** `MD| Pressure` additionally contains
the kinetic term and must not be used for a training set — it would bake kinetic
energy into the potential and give systematically wrong pressures at other
temperatures. For the reference run, frame 1:

$$P_\text{total} = \frac{2}{3}\frac{E_\text{kin}}{V} + \frac{1}{3}\mathrm{Tr}\,\sigma
= 14134 + 2345 = 16479\ \text{bar}$$

against `MD| Pressure = 16480.9 bar` — the 1.6 bar residual is rounding in the
printed values.

## Dropped frames

Three kinds of frame are discarded, and the counts are always reported:

| Reason | Why |
|---|---|
| SCF not converged | the forces are garbage; a gap is better than bad labels |
| incomplete block | a job killed mid-step leaves a truncated frame |
| composition changed | guards against block misalignment, see below |

For VASP a frame is also dropped when it has **no cell block of its own**. The
cell is never inherited from the previous frame: in a fixed-cell run the two are
identical so the bug would be invisible, and it would then produce silently
wrong data the first time someone ran a variable-cell job.

The composition check is not really about the system changing. The two xyz
blocks CP2K prints per step — coordinates then forces — are *byte-for-byte
indistinguishable*: same atom count, same `i = …, time = …, E = …` comment
line, same element column. Only their order tells them apart. If some warning
is printed inside a block, the element column stops holding element symbols and
the composition check catches it.

If a large fraction of frames is dropped, that is a signal about the run, not
about ferro — a 5000-frame trajectory reduced to 2000 usually means the SCF
settings need attention.

A `WARNING: … frame(s) print their blocks at a different offset than the first`
means extra output is interleaved between the blocks. Nothing was necessarily
lost — the scan is offset-independent — but it is worth checking a few frames by
hand.

## What is not done here

- **No filtering.** Removing frames by force/stress magnitude or by geometry is
  [`ferro dataset filter`](filter.md).
- **No cross-run merging or shuffling.** Combining datasets from *different*
  runs and resizing sets is [`ferro dataset merge`](merge.md). The files of one
  directory are a different case — they are one run, and `collect` reassembles
  them.
- **No extxyz / NEP output.** The pipeline's intermediate format is the DeePMD
  directory; GPUMD's `train.xyz` is an export at the end of the chain, not at
  the start.
- **None of these formats is registered with `ferro convert`.** `.out` is far
  too generic to be claimed for CP2K, and `OUTCAR` has no extension at all, so
  the AIMD readers are reachable only through `ferro dataset collect`.
- **No ML force-field OUTCARs.** A VASP run driven by its machine-learned force
  field prints `free  energy ML TOTEN` and `ML FORCE` instead, and its block
  layout differs by more than the names. Without a sample to check against,
  guessing would be worse than declining.
