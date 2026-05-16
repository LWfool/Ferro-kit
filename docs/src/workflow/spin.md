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

\\[ n_\text{unpaired} = \operatorname{round}\!\left( \left| \sum_i m_i \right| \right), \qquad 2S+1 = n_\text{unpaired} + 1 \\]

This reflects an actual (DFT- or user-supplied) spin state and is preferred whenever magmom data is present.

## Strategy 2 — Oxidation state + Hund's rule

For ionic solids with no magmom data, `assign_oxidation_states` assigns formal oxidation states by electronegativity rules and charge balance:

1. The most electronegative element that has a negative common oxidation state is the **anion**; it takes its most-negative common state (O → −2, F → −1, S → −2, …).
2. Remaining elements are **cations**.  The combination of their common positive oxidation states that satisfies overall charge balance is selected (highest-oxidation solution preferred when ambiguous).

Each ion's unpaired count is then:

- **Transition metals** — d-electron count \\( n_d = \text{group} - \text{oxidation state} \\); high-spin filling of 5 d-orbitals (Hund's rule):
  \\[ n_\text{unpaired} = \begin{cases} n_d & n_d \le 5 \\\\ 10 - n_d & n_d > 5 \end{cases} \\]
- **Main group** — valence electrons after ionization filled into the s/p shell; closed-shell ions give 0.

Contributions are summed over all sites (ferromagnetic assumption → upper bound).

## Strategy 3 — Electron-count parity

With no usable structural information, only a bound is given from the total electron count \\( N_e = \sum_i Z_i - q \\):

- \\( N_e \\) odd → at least one unpaired electron (doublet, 2S+1 = 2)
- \\( N_e \\) even → singlet assumed (2S+1 = 1)

This is always applied as a final sanity check: if the multiplicity from strategy 1 or 2 has the wrong parity relative to \\( N_e \\), a warning is emitted and the parity bound is used instead.

---

## Worked examples

### ZnP₂O₆ — diamagnetic

| Ion | Configuration | Unpaired |
|---|---|---|
| Zn²⁺ | group 12, \\(n_d = 12-2 = 10\\) → d¹⁰ | 0 |
| P⁵⁺ | main group, 5 − 5 = 0 valence e⁻ | 0 |
| O²⁻ | 6 + 2 = 8 → closed octet | 0 |

Oxidation states resolve uniquely (Zn fixed at +2 ⇒ P = +5).  Total = **0 unpaired → multiplicity 1**.

### MnS — high-spin d⁵

S (electronegativity 2.58) is the anion at −2 ⇒ Mn = +2.
Mn²⁺: group 7, \\(n_d = 7 - 2 = 5\\) → d⁵ high-spin → **5 unpaired → multiplicity 6**.

### Fe₂O₃

O₃ = −6 ⇒ 2 Fe = +6 ⇒ Fe = +3.  Fe³⁺: \\(n_d = 8 - 3 = 5\\) → 5 unpaired each → **10 total → multiplicity 11**.

---

## Usage

```bash
# Auto-guess for a transition-metal oxide (CP2K)
fe-job -s cp2k -i Fe2O3.cif --auto-spin --smear

# Auto-guess for QE
fe-job -s qe -i Fe2O3.cif --auto-spin --kpoints 4 4 4

# Manual override (highest priority — disables auto-spin)
fe-job -s gaussian -i radical.xyz --charge 0 --multiplicity 2
```

Priority: explicit `--multiplicity` > `--auto-spin` (or builder default) > value from the input file.  In the CP2K and QE builders `auto_spin` is on by default; passing `--multiplicity` disables it so the manual value is respected.

---

## Limitations

These are inherent to formal-charge / parity methods and are reported as warnings:

- **Transition metals are always estimated high-spin.**  Low-spin complexes (strong-field ligands) and the geometry-dependent crossover are not detected — verify with DFT.
- **Multi-centre magnetic coupling** (ferro- vs antiferromagnetic) cannot be inferred from structure; the sum is an upper bound.
- **Purely covalent molecules** fall back to the parity bound.  O₂, for example, is reported as a singlet — its triplet ground state arises from π* orbital degeneracy, which requires molecular-orbital theory.
- **Lanthanide f-electrons** are not counted.
- **Covalent transition-metal complexes** (organometallics) do not satisfy the ionic assumption.

The estimate is an *initial guess*, not a substitute for an electronic-structure calculation.

---

## Related

- [Job Builders](job-builders.md) — Gaussian / CP2K / QE input generation
- [CLI Reference: `fe-job`](../cli-reference.md#fe-job)
