# Averaged Charge-Density SDF (`chg-sdf`)

`fe-cube -m chg_sdf` computes the **averaged electron-density spatial distribution function** over a set of Qn-type network clusters.  Given a collection of QE `pp.x` charge-density cube files (one per MD snapshot), the tool:

1. Identifies every Qn cluster in each snapshot using the same Kabsch-alignment pipeline as `fe-cube -m sdf`.
2. Extracts a cubic sub-grid of the electron density centred on the cluster anchor atom, interpolating the raw charge-density grid with trilinear interpolation (PBC-safe).
3. Rotates the sub-grid to align with the first-encountered reference cluster via the Kabsch rotation matrix, using "pull" interpolation (`output[r] = input[R⁻¹ r]`).
4. Accumulates rotated sub-grids; divides by the cluster count to obtain the averaged density.
5. Writes one Gaussian cube file per signature family.

The result reveals the orientationally averaged electron-density distribution around a specific local structural motif (e.g. Q2 PO₄ unit in a phosphate glass).

---

## Prerequisites

Charge-density cube files must be produced by QE `pp.x` (or any program that outputs the Gaussian cube format with atom coordinates in Bohr and density in `e/Bohr³`).  One cube file per MD frame is expected.

---

## Usage

```bash
# Basic: average charge density around Q2 clusters (P–O system)
fe-cube -m chg_sdf \
    --cubes frame_000.cube frame_001.cube frame_002.cube \
    --qn 2 --former P --ligand O --cutoff-fl 2.4 \
    -o chg_sdf

# With modifier cation (e.g. Zn2+ in Zn-phosphate glass)
fe-cube -m chg_sdf \
    --cubes *.cube \
    --qn 2 --former P --ligand O --cutoff-fl 2.4 \
    --modifier Zn --cutoff-ml 2.8 \
    --chg-padding 6.0 \
    -o chg_sdf_Q2

# Tighter sub-grid / multiple threads
fe-cube -m chg_sdf \
    --cubes *.cube \
    --qn 3 --former P --ligand O \
    --chg-padding 5.0 --rmsd-warn 0.3 \
    --ncore 8 \
    -o q3_density
```

---

## Parameters

| Flag | Default | Description |
|---|---|---|
| `--cubes <files…>` | (required) | One or more QE pp.x cube files; each file = one MD frame |
| `--qn` | `3` | Target Qn cluster level (0, 1, 2, or 3) |
| `--former` | `P` | Network-former element symbol |
| `--ligand` | `O` | Bridging-ligand element symbol |
| `--cutoff-fl` | `2.4` | Former–ligand bond cutoff [Å] |
| `--modifier` | (none) | Modifier cation element; omit if not needed |
| `--cutoff-ml` | `2.8` | Modifier–ligand cutoff [Å] |
| `--chg-padding` | `6.0` | Sub-grid margin around the cluster anchor [Å] |
| `--rmsd-warn` | `0.5` | RMSD threshold for alignment quality warnings [Å] |
| `--ncore` | all | Parallel threads |
| `-o` | `chg_sdf` | Output file stem |

> **Note**: `--chg-padding` controls the sub-grid size, not the atom SDF grid.  Larger values capture more surrounding density but increase memory and runtime.

---

## Output

For a single signature family:
```
chg_sdf_Q2.cube
```

For multiple families (distinct cluster topologies at the same Qn):
```
chg_sdf_fam0_Q2.cube
chg_sdf_fam1_Q2.cube
...
```

The cube file contains:
- **Atom positions**: the reference cluster atoms in a local Cartesian frame (anchor at origin).
- **Density values**: averaged `ρ_phys × V_cell` in internal units.  To convert to physical electron density [e/Å³], divide by the voxel volume in Å³.

Console output reports per-family cluster count and alignment RMSD statistics:

```
Family "P-O-O-P"  (142 clusters, RMSD mean=0.031 max=0.187 Å, 0 warnings) → chg_sdf_Q2.cube
ChgSDF Q2 done: 50 frames, 142 clusters total, 1 cube file written
```

---

## Algorithm Details

### Cluster identification
Identical to `fe-cube -m sdf`: finds all Qn-type connected components of former atoms sharing bridging ligands, using Union-Find on the former–ligand–former connectivity graph.

### Sub-grid extraction
For each cluster anchor position `c` (Cartesian, Å):

```
sub[ix, iy, iz] = ρ_charge( c + Δr )
```

where `Δr = (ix − half_n, iy − half_n, iz − half_n) × δ_voxel` and `ρ_charge` is evaluated via trilinear interpolation with PBC folding on the original charge-density grid.  The sub-grid has shape `n × n × n` with `n = 2 × ⌈padding / δ_voxel⌉ + 1`.

### Grid rotation
For each subsequent cluster the Kabsch rotation matrix `R` is computed from the atom positions (same as the atom SDF).  The sub-grid is rotated using "pull" interpolation:

```
output[r_out] = input[ R⁻¹ · (r_out − centre) + centre ]
```

`R⁻¹ = Rᵀ` for a proper rotation.  Points outside the sub-grid boundary return 0.

### Accumulation
```
ρ_avg = (1 / N) Σᵢ ρᵢ_rotated
```

---

## Related Commands

- [`fe-cube -m sdf`](cube-sdf.md) — atom-type SDF (no charge density required)
- [`fe-bader`](../cli-reference.md#fe-bader) — Bader charge decomposition from CHGCAR / cube files
