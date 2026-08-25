# Collecting AIMD Output

`ferro dataset collect` turns *ab initio* MD output into
[DeePMD-kit](https://docs.deepmodeling.com/projects/deepmd/) system
directories — the starting point for training a machine-learning potential.

```bash
ferro dataset collect -i total.out                 # -> ./total/
ferro dataset collect -i run*/total.out -o data    # -> data/run1_total/, ...
```

Currently only **CP2K MD output** is read. VASP and Quantum ESPRESSO are
planned; the parsing layer and the writing layer are kept separate so adding a
source means writing only the former.

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

## Output layout

One system directory per input file:

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

`<name>` is the input file stem. CP2K logs are routinely all called
`total.out`, so when several inputs share a stem the parent directory is
prefixed — `run1_total`, `run2_total` — instead of one silently overwriting the
next.

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
  `ferro dataset filter` (not implemented yet).
- **No merging or shuffling.** Combining same-composition datasets and resizing
  sets is `ferro dataset merge` (not implemented yet).
- **No extxyz / NEP output.** The pipeline's intermediate format is the DeePMD
  directory; GPUMD's `train.xyz` is an export at the end of the chain, not at
  the start.
- **`.out` is not registered with `ferro convert`.** That extension is far too
  generic to be claimed for CP2K, so the CP2K MD reader is reachable only
  through `ferro dataset collect`.
