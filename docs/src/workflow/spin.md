# Unpaired-Electron / Spin-Multiplicity Estimation

When generating a QC input file you must specify the system's **spin multiplicity** (2S+1).  Guessing wrong produces a calculation that either fails to converge or converges to the wrong electronic state.  `ferro` can estimate the multiplicity directly from the structure via `ferro_core::guess_spin`, exposed on every job builder through the `--auto-spin` flag.

The estimator uses three strategies in decreasing order of reliability and falls back automatically:

1. **Magnetic-moment sum** — most reliable
2. **Oxidation state + Hund's rule** — for ionic solids
3. **Electron-count parity** — universal lower bound

The result is always reconciled against the total electron count; an inconsistency triggers a warning and a fall back to the parity bound.

---

## Strategy 1 — Magnetic-moment sum

If any atom carries a `magmom` (read from a QE input, extxyz, or VASP OUTCAR), the number of unpaired electrons is

$$ n_\text{unpaired} = \operatorname{round}\!\left( \left| \sum_i m_i \right| \right), \qquad 2S+1 = n_\text{unpaired} + 1 $$

This reflects an actual (DFT- or user-supplied) spin state and is preferred whenever magmom data is present.

## Strategy 2 — Oxidation state + Hund's rule

For ionic solids with no magmom data, `assign_oxidation_states` assigns formal oxidation states by electronegativity rules and charge balance:

1. The most electronegative element that has a negative common oxidation state is the **anion**; it takes its most-negative common state (O → −2, F → −1, S → −2, …).
2. Remaining elements are **cations**.  The combination of their common positive oxidation states that satisfies overall charge balance is selected (highest-oxidation solution preferred when ambiguous).

Each ion's unpaired count is then:

- **3d transition metals** (Sc–Zn) — d-electron count $ n_d = \text{group} - \text{oxidation state} $; high-spin filling of 5 d-orbitals (Hund's rule):
  $$ n_\text{unpaired} = \begin{cases} n_d & n_d \le 5 \\\\ 10 - n_d & n_d > 5 \end{cases} $$
- **4d / 5d transition metals** (Y–Cd, Hf–Hg) — same $ n_d $, **low-spin**: the three $t_{2g}$ orbitals fill before $e_g$, giving 0, 1, 2, 3, 2, 1, 0, 1, 0, 1, 0 for $d^0$…$d^{10}$. $d^8$ is taken as square-planar (0 unpaired): $Pd^{2+}$, $Pt^{2+}$ and $Au^{3+}$ are almost always square-planar. Their ligand-field splitting is large enough that high-spin 4d/5d ions are rare.
- **Lanthanides** (La–Lu) — the ion keeps its electrons beyond the Xe core in 4f, $ n_f = Z - 54 - \text{oxidation state} $, filled into 7 f-orbitals by Hund's rule (free-ion high spin). $Ce^{4+}$ $f^0$ → 0, $Ce^{3+}$ $f^1$ → 1, $Gd^{3+}$ and $Eu^{2+}$ $f^7$ → 7, $Eu^{3+}$ $f^6$ → 6. Spin–orbit coupling does not enter: $Eu^{3+}$ has a $J = 0$ ground state but $S = 3$, and a scalar-relativistic calculation needs the multiplicity $2S+1 = 7$.
- **Main group** — valence electrons after ionization filled into the s/p shell; closed-shell ions give 0.

Contributions are summed over all sites (ferromagnetic assumption → upper bound).

## Strategy 3 — Electron-count parity

With no usable structural information, only a bound is given from the total electron count $ N_e = \sum_i Z_i - q $:

- $ N_e $ odd → at least one unpaired electron (doublet, 2S+1 = 2)
- $ N_e $ even → singlet assumed (2S+1 = 1)

This is always applied as a final sanity check: if the multiplicity from strategy 1 or 2 has the wrong parity relative to $ N_e $, a warning is emitted and the parity bound is used instead.

---

## Worked examples

### $ZnP_2O_6$ — diamagnetic

| Ion | Configuration | Unpaired |
|---|---|---|
| $Zn^{2+}$ | group 12, $n_d = 12-2 = 10$ → $d^{10}$ | 0 |
| $P^{5+}$ | main group, 5 − 5 = 0 valence $e^-$ | 0 |
| $O^{2-}$ | 6 + 2 = 8 → closed octet | 0 |

Oxidation states resolve uniquely (Zn fixed at +2 ⇒ P = +5).  Total = **0 unpaired → multiplicity 1**.

### MnS — high-spin $d^5$

S (electronegativity 2.58) is the anion at −2 ⇒ Mn = +2.
$Mn^{2+}$: group 7, $n_d = 7 - 2 = 5$ → $d^5$ high-spin → **5 unpaired → multiplicity 6**.

### $Fe_2O_3$

$O_3$ = −6 ⇒ 2 Fe = +6 ⇒ Fe = +3.  $Fe^{3+}$: $n_d = 8 - 3 = 5$ → 5 unpaired each → **10 total → multiplicity 11**.

### $CeO_2$ and $Ce_2O_3$ — 4f

$CeO_2$: Ce = +4, $n_f = 58 - 54 - 4 = 0$ → **multiplicity 1**.
$Ce_2O_3$: Ce = +3, $n_f = 1$ → 1 unpaired each → **2 total → multiplicity 3**.

---

## Usage

```bash
# Auto-guess for a transition-metal oxide (CP2K)
ferro job -s cp2k -i Fe2O3.cif --auto-spin --smear

# Auto-guess for QE
ferro job -s qe -i Fe2O3.cif --auto-spin --kpoints 4 4 4

# Manual override (highest priority — disables auto-spin)
ferro job -s gaussian -i radical.xyz --charge 0 --multiplicity 2
```

Priority: explicit `--multiplicity` > `--auto-spin` (or builder default) > value from the input file.  In the CP2K and QE builders `auto_spin` is on by default; passing `--multiplicity` disables it so the manual value is respected.

---

## Limitations

These are inherent to formal-charge / parity methods and are reported as warnings:

- **Spin state is fixed by row, not by ligand field.**  3d ions are always high-spin and 4d/5d ions always low-spin. Low-spin 3d complexes (strong-field ligands, e.g. $Co^{3+}$ $d^6$), the geometry-dependent crossover and octahedral 4d/5d $d^8$ are not detected — verify with DFT.
- **Multi-centre magnetic coupling** (ferro- vs antiferromagnetic) cannot be inferred from structure; the sum is an upper bound.
- **Purely covalent molecules** fall back to the parity bound.  $O_2$, for example, is reported as a singlet — its triplet ground state arises from π* orbital degeneracy, which requires molecular-orbital theory.
- **Lanthanides assume 4f in the valence.**  Pseudopotentials or ECPs that freeze 4f in the core (VASP `Ln_3`, Stuttgart large-core ECPs) need a different count, and the electron parity changes with them; give `--multiplicity` by hand. All CP2K GTH lanthanide potentials that `ferro job` selects keep 4f in the valence. Actinides are not covered.
- **One integer oxidation state per element.**  Mixed-valence compounds such as $Fe_3O_4$ ($Fe^{2+}$ + 2 $Fe^{3+}$) cannot be balanced; the estimator says so and falls back to the parity bound.
- **Covalent transition-metal complexes** (organometallics) do not satisfy the ionic assumption.

The estimate is an *initial guess*, not a substitute for an electronic-structure calculation.

---

## Related

- [Job Builders](job-builders.md) — Gaussian / CP2K / QE input generation
- [CLI Reference: `ferro job`](../cli-reference.md#ferro job)
