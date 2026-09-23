# Glass Network Analysis (`ferro net`)

`ferro net` builds the glass network topology frame by frame over an MD trajectory, assigns every
atom a **structured type**, and accumulates the distributions over the whole trajectory.  The output is six long tables plus an optional labelled trajectory.

| Quantity | Table | Meaning |
|---|---|---|
| **Structural composition** | `composition` | one species per row: `P-Q2` `Al_4` `O_b` `Zn_4`, each as a fraction of its own element |
| Qn distribution | `qn` | the **homonuclear** connection count $n$ of a Qn network former (P–O–P), the $n$ of the literature's $Q^n_m$ |
| Heteronuclear bridge decomposition | `qn_partner` | the table above split one dimension further by partner element, i.e. $Q^n(m\mathrm{Al})$ |
| Ligand classification | `ligand_type` | `O_f` / `O_n` / `O_b` / `O_t`, with the partner elements given as data columns |
| Coordination number | `coordination` | the total coordination-number distribution of network formers and modifiers |
| Linkage statistics | `linkage` | the ligand element of every bridge and the site state at both of its ends |

`composition` is a **summary of the other tables**, not a new measurement: read it to take in the structural
composition of the glass at a glance, and go to the source tables when the partner decomposition or the state at both ends is needed.

> **Qn is reported only for Qn elements.** Qn is the notation for tetrahedral network formers, so by
> default only `B` / `P` / `Si` appear in the `qn` and `qn_partner` tables; **network formers such as Al
> are characterised by coordination number** (written Al[4] / Al[5] / Al[6] in the literature).  They still
> take part in the bridging-oxygen decision, in `ligand_type` and in `linkage` — what they leave is only the *rows* of those two tables.  `--qn` replaces this list as a whole.

> **Command shape.** `ferro net qn` and `ferro net type` have been merged into the leaf command
> `ferro net`; the export is now just the `--export-traj` switch.  The old labels (`P0` / `Of` / `On_P` /
> `Ob_P_P` / `X` / `Zn_f`) have all been replaced, with no read-compatibility layer.

---

## Theory

### Deciding which atoms are bonded

A cutoff radius is the bonding criterion: two atoms are bonded when the distance between them is below
the cutoff.  Distances use the minimum-image convention, which supports orthorhombic and triclinic cells:

$$d_{ij} = \bigl|\mathbf{r}_{ij} - \mathbf{M} \cdot \text{round}\!\left(\mathbf{M}^{-1}\mathbf{r}_{ij}\right)\bigr|$$

The analysis requires a cell (PBC) on every frame.  Frames without one are skipped; when no frame is left, that input errors out and is skipped.

### The three roles

The parameters split the elements into three classes, and everything below rests on that split:

| Role | Where it comes from | Takes part in |
|---|---|---|
| **network former** | the left-hand side of `--<F>-<L>=<Å>` | number of bridges, coordination number, ligand classification, linkage statistics |
| **ligand** | the right-hand side of `--<F>-<L>=<Å>` | ligand classification |
| **modifier** | declared with `--modifier`; its cutoffs also use `--<M>-<L>=<Å>` | coordination number **only** |

Keeping modifiers out of the ligand classification is the entire reason this parameter exists: a Zn sitting
next to a non-bridging oxygen must not turn that oxygen into a bridging one.  Declaring a modifier without
giving it a cutoff errors out at once — letting it pass silently would make it count as a network former, and the whole oxygen classification would shift with it.

### Ligand classification

For every ligand atom $k$, count its network-former neighbours $n_k$:

| $n_k$ | Label | Meaning |
|---|---|---|
| 0 | `O_f` | Free — bonded to no network former |
| 1 | `O_n` | Non-bridging ligand |
| 2 | `O_b` | Bridging ligand |
| ≥3 | `O_t` | Tricluster — a three-coordinate ligand, common in Al-bearing glasses |

**The label carries no partner element.** Both `P–O–P` and `P–O–Al` are labelled `O_b`; the two are told
apart by the `former_a` / `former_b` **data columns** of the `oxy` table.  The label stays plain and the
statistics stay complete, which decouples the two: the label is for structure files, the partner decomposition is for analysis.

### Two quantities that must be kept apart: number of bridges and coordination number

For every network-former atom $i$:

$$n_i^\text{bridge} = \bigl|\{k \in \text{neighbors}(i) : n_k \ge 2\}\bigr|,
\qquad
\mathrm{CN}_i = \bigl|\text{neighbors}(i)\bigr|$$

**The number of bridges counts bridging ligands only; the coordination number counts every ligand inside
the cutoff**, non-bridging ones included.  A three-coordinate ligand does count towards the number of
bridges (it does connect that former into the network), at the price of $\sum(\text{bridges}) \neq 2 \times |\text{O\_b}|$ — one `O_t` is counted once by each of its three sides.

> **Their being equal is a property of the system, not a definition.** A network former with no
> non-bridging ligand makes $n^\text{bridge} = \mathrm{CN}$ hold everywhere.  That is exactly the case for
> Al in the reference trajectory (the `ligand_type` table has no `O_n, Al` row at all), so the two
> distribution tables agree bin by bin.  Switch to a system where Al carries non-bridging oxygens and the two part ways at once.  Whenever you read "the coordination number of Al", take it from `coordination`.

### Qn is only for Qn network formers

$Q^n$ is the notation for tetrahedral network formers: on a site whose coordination number is essentially
fixed, one number says everything about how it is connected.  Al does not meet that premise — its coordination number is itself one of the quantities to be reported.  Hence:

| | Appears in `qn` / `qn_partner` | Number in the label | Characterised by |
|---|---|---|---|
| Qn network former (B, P, Si by default) | yes | $n$ (homonuclear connection count) | the Qn distribution |
| other network formers (Al, …) | **no** | **coordination number** | the `coordination` table |
| modifiers (Zn, …) | no | no suffix | the `coordination` table |

The list is a **default**: `--qn Si,Al` replaces it as a whole rather than adding to it — whether an element
counts as a Qn network former is a property of the system: in aluminosilicates people do report $Q^n(m\mathrm{Si})$ for Al.

### $n$ counts homonuclear bridges only — the statement used in the literature

This is **the column most easily misread**; please finish this section before using `qn`.

In the extended notation $Q^n_m$ of the literature (written $Q^n(m\mathrm{Al})$ for aluminophosphates,
$Q^n(m\mathrm{B})$ for borophosphates, and $P^n_{m\mathrm{Al},x\mathrm{B}}$ with several heteronuclear partners):

| Symbol | What it counts |
|---|---|
| $n$ | **homonuclear** bridging connections: P–O–P |
| $m_X$ | **heteronuclear** bridging connections: P–O–Al, P–O–B … |
| total bridging connections | $n + \sum_X m_X$ — **not** $n$ |

> Verbatim from arXiv 2510.13545: *"'n' denotes the total number of bridging oxygen
> connections involving P–O–P and P–O–Si bonds, and 'm' and 'x' specify the number
> of aluminum and boron connections with phosphate."*
>
> The aluminophosphate literature points out explicitly that reading "$n$ is the total number of bridging
> oxygens" is a **known misconception** in this field.  It arises because $n$ and $m$ come from two
> different NMR experiments (homonuclear J / double-quantum vs heteronuclear REAPDOR / TRAPDOR /
> HETCOR) and are two independent counts by nature.  Species with $m>n$ such as $Q^0(3\mathrm{B})$ and $Q^1(2\mathrm{Al})$ are real and common (the former is 29% at high boron content), and cannot be explained under the reading that "$m$ is a subset of $n$".

`ferro net` has followed this statement since 2026-08.  Its effect on the reference trajectory `43Z43P15A`:

| | Old statement (total bridges) | Literature statement |
|---|---|---|
| `mean_qn` | 2.40 | **0.95** (with $m_\mathrm{Al}$=1.45 alongside) |
| `P-Q3` fraction | 40.4% | **0.27%** |

A factor of 130.  Under the old statement a reader would take "`P-Q3` is 40%" to mean a highly
cross-linked phosphate network, while the P in this glass carries on average only 0.95 P–O–P bridges — a structure dominated by dimers and short chains.

**The total bridging-oxygen count is still available**: the `[inputs]` block gives `mean_n_bo` next to `mean_qn`.

- `qn=1, m_Al=2` → $Q^1(2\mathrm{Al})$, three bridges in total
- `qn` is the **marginal** of `qn_partner` over the `m_` columns (count and fraction close exactly)

**The partner decomposition is not encoded into `label`.** With several network formers the label would
grow into `P-Q1(2Al,1B)`, whereas the `qn` table exists precisely to be scanned at a glance; the `m_<X>`
columns are what filtering and plotting are for.  The literature itself uses subscript positions in $P^n_{m\mathrm{Al},x\mathrm{B}}$ rather than a string of parentheses.

**There is no `m_` column for the former's own element**: under the new statement `m_P` is identically
equal to `qn`, a duplicate column.  It is fully determined by `qn` and is not an independent grouping dimension, so dropping it leaves the closure above intact.

Note that Al is a **partner here, not a subject**: it has no row of its own, but the `m_Al` column is there as usual.

The two are **two tables rather than one plus a `groupby`**: the plain Qn distribution is the primary output
and must be readable the moment the file is opened, without an aggregation first.  This is the same
criterion that made `average` a table of its own — a different granularity deserves a separate table.  `sd`
likewise has to be accumulated separately rather than summed: the variance of a sum of correlated terms is
not the sum of the variances.  Measured on the reference trajectory, the `qn=1` row of P has a correct `sd`
of **0** (the P–O–P backbone does not change from frame to frame), while summing the `sd` of the four corresponding `qn_partner` rows gives **7.598e-3** — summing would forge a visible fluctuation for a quantity that does not move at all.  The components moving while their sum stays put is exactly what cancellation between correlated terms looks like.

**A tricluster ligand is counted by connections.** A ligand bonded to three network formers connects the
site in question to **two** partners, so it contributes 2 rather than 1 — the literature counts *connections*,
not bridging oxygens.  As a result $n + \sum_X m_X$ can **exceed** the bridging-oxygen count (`mean_n_bo`),
and the excess is exactly the tricluster bridges.  This is also what lets $Q^n_m$ work as an unambiguous label: counting by bridging oxygens would leave $\sum m$ short.

This decomposition **cannot be recovered from the `linkage` table**: two P–O–Al bridges may come from one
P with $m_\mathrm{Al}=2$, or from two P with $m_\mathrm{Al}=1$.  `linkage` counts bridges, `qn` counts atoms.

### Linkage statistics

Every ligand bonded to ≥2 network formers records one entry:

$$(\underbrace{\text{element},\; n_\text{bridge},\; \text{CN}}_{\text{A end}}),\;
  (\underbrace{\text{element},\; n_\text{bridge},\; \text{CN}}_{\text{B end}}),\;
  \text{ligand element},\; n_\text{formers}$$

**Both** ends carry all three fields.  This differs from common implementations, which store a single
number per element (Qn for P, CN for Al) and therefore cannot answer "how many bridges does a four-coordinate Al have".

- **The `linkage` display column**: a human-readable form such as `Al_4-O-P_2`, where the number follows
  each one's own convention (Qn for Qn network formers, coordination number for the rest), in the same
  vocabulary as the labels of the exported trajectory.  **Filter on the numeric columns**; do not take this one apart with a regex.
- **The `ligand` column**: the element of the atom in the middle of the bridge.  With several ligand species
  `Al-O-P` and `Al-F-P` are two different bridges and never share a row — merging their counts would report a number no experiment can correspond to.
- **Canonical half**: a linkage has no direction, so the two ends are sorted by `(element, homonuclear connection count, CN)` with the smaller first, and each pair is stored once.
  A **row sum is therefore not that site's total participation**; to get participation, count both the `_a` and the `_b` columns.
- **The `n_formers` column**: 2 for an ordinary bridging oxygen; a three-coordinate ligand expands into $C(3,2)=3$ rows, each marked 3.
  They are not discarded, because "what a three-coordinate oxygen connects to" is exactly what one wants to study in Al-bearing systems.
- **`qn_a/b` are homonuclear connection counts**, defined for every network former — non-Qn ones included
  (the `qn_a` of Al is its Al–O–Al count).  They have the **same origin** as the number in the `linkage`
  display column: the 2 of `P_2` is read from `qn_a`, the 4 of `Al_4` from `cn_a`, so the label and the numeric columns can never contradict each other.

### Statistical statement: `fraction` and `sd`

- `count`: the sum of the counts over the whole trajectory (integer, exact)
- `fraction`: the **mean of the per-frame fractions**.  A frame without that bin counts as 0, not as a missing value
- `sd`: the **sample standard deviation** of the same series (ddof = 1), empty for a single frame

> **`sd` is not a standard error.** Neighbouring MD frames are strongly correlated, so it neither shrinks as
> $1/\sqrt{N}$ nor estimates the physical fluctuation.  Read it as "how much this number swings between
> snapshots".  Real error bars call for block averaging, which this analysis does not provide.

The implementation uses Welford's recurrence rather than $\Sigma f$ / $\Sigma f^2$: the latter produces
~2e-10 of cancellation noise on a constant series, and a spurious wobble in the `sd` column would be read as physics.

---

## Quick start

```bash
# P2O5 glass, P-O cutoff 2.4 Å
ferro net -i traj.lammpstrj --P-O=2.4

# aluminophosphate + Zn modifier
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn

# batch: several trajectories in one run, stacked into one csv with a file column
ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan

# last 500 frames only, 8 threads
ferro net -i traj.lammpstrj --P-O=2.4 --last-n 500 --ncore 8

# also export the labelled trajectory
ferro net -i traj.lammpstrj --P-O=2.4 --export-traj
ferro net -i traj.lammpstrj --P-O=2.4 --export-traj extxyz
```

---

## Command-line parameters

### Required

At least one pair parameter of the form `--<Former>-<Ligand>=<cutoff>` (element symbols capitalised):

```
--P-O=2.4        cutoff 2.4 Å between P (network former) and O (ligand)
--Si-O=1.8       Si-O cutoff 1.8 Å
--Al-O=2.4       Al-O cutoff 2.4 Å
--Al-F=2.1       one network former may have several ligand species
```

The element pair lives in the **parameter name**, which clap cannot model, so `main` strips them out of argv before parsing.

### Optional

| Parameter | Default | Notes |
|---|---|---|
| `-i <FILE>...` | — | input trajectories; takes several values and expands globs itself (quote them).  Prints the help when omitted |
| `-o <DIR>` | current directory | the six tables and the `--export-traj` trajectory all go here.  Asks first when the directory does not exist; `--mkdir` skips the question |
| `-s <SUFFIX>` | — | batch marker: `network_<table>_<suffix>.csv` |
| `--mkdir` | off | create the `-o` directory without asking (required in non-interactive environments) |
| `--last-n N` | all frames | use only the last N frames |
| `--ncore N` | all cores | number of parallel threads |
| `--metal-units` | off | LAMMPS metal units.  The statistics read neither velocities nor forces, so this **only affects `--export-traj extxyz`** |
| `--modifier E,E` | — | comma-separated elements that count towards coordination number only.  Their cutoffs must be given as well |
| `--qn E,E` | `B,P,Si` | comma-separated network formers to report Qn for.  **Replaces** the default list rather than adding to it |
| `--export-traj [FMT]` | — | also write a labelled trajectory, `lammpstrj` (default) or `extxyz` |

---

## Output

Six CSVs, each with a `file` column.  The `#` comment block of every file holds the shared parameters, the
`[inputs]` list **and that table's own column-by-column description** (`pandas.read_csv(comment="#")` drops
the whole block automatically).  The table below therefore only says which file holds what; the meaning of the columns is in the file itself.

| File | One row per | Columns |
|---|---|---|
| `network_composition.csv` | species | `label, element, count, fraction, sd` |
| `network_qn.csv` | (network former, Qn) | `label, former, qn, count, fraction, sd` |
| `network_qn_partner.csv` | (network former, Qn, partner decomposition) | `label, former, qn, m_<X>…, count, fraction, sd` |
| `network_ligand_type.csv` | (ligand type, partner pair) | `label, type, former_a, former_b, count, fraction, sd` |
| `network_coordination.csv` | (element, coordination number) | `element, cn, count, fraction, sd` |
| `network_linkage.csv` | (ligand, state at both ends) | `linkage, ligand, elem_a, qn_a, cn_a, elem_b, qn_b, cn_b, n_formers, count, fraction, sd` |

The `label` column is the human-readable anchor; the numeric columns `former` / `qn` / `cn` are what
filtering and plotting use — both are kept rather than one or the other, otherwise "select Qn ≥ 3" would mean slicing strings.

The value columns of the distribution tables mean different things (Qn / type label / coordination number),
so merging them would make `groupby` meaningless and each gets its own table.  `qn` vs `qn_partner`, and `average`, are all cases of **a different granularity deserving a separate table**.

**When there is no Qn network former the first two files are not written at all**, and the reason is printed
to the screen.  A CSV with nothing but a header reads as "measured, and the result was zero", while in fact nothing was measured.

In batch mode, inputs with different element sets take the **union of the columns; what is missing stays
empty (NaN), never padded with zeros, never interpolated** — a system without Al has an empty `m_Al` column, not 0.

### Example: `network_composition.csv`

```csv
file,label,element,count,fraction,sd
43Z43P15A,P-Q0,P,465,2.500000e-1,0.000000e0
43Z43P15A,P-Q1,P,1035,5.564516e-1,0.000000e0
43Z43P15A,P-Q2,P,355,1.908602e-1,0.000000e0
43Z43P15A,P-Q3,P,5,2.688172e-3,0.000000e0
43Z43P15A,Al_4,Al,590,8.939394e-1,1.417294e-2
43Z43P15A,Al_5,Al,63,9.545455e-2,1.570943e-2
43Z43P15A,Al_6,Al,7,1.060606e-2,4.149413e-3
43Z43P15A,Zn_3,Zn,87,9.354839e-2,4.906927e-2
43Z43P15A,Zn_4,Zn,675,7.258065e-1,2.434243e-2
43Z43P15A,Zn_5,Zn,159,1.709677e-1,3.921416e-2
43Z43P15A,Zn_6,Zn,8,8.602151e-3,2.944745e-3
43Z43P15A,Zn_7,Zn,1,1.075269e-3,2.404374e-3
43Z43P15A,O_n,O,2985,4.543379e-1,1.203302e-3
43Z43P15A,O_b,O,3583,5.453577e-1,1.154167e-3
43Z43P15A,O_t,O,2,3.044140e-4,4.168360e-4
```

> The `sd` of 0 on the four P rows is not a defect: the topology of the P–O–P backbone does not change
> across these 5 frames.  All the fluctuation sits on the P–O–Al and coordination-number side (look at the
> `sd` of `Al_*` and `Zn_*`).  The old statement mixed the rigid backbone and the fluctuation into one number, which is why these rows used to have a non-zero `sd`.

**The denominator is always the atom count of that element**: Q2 out of all P, `Al_4` out of all Al, `O_b`
out of all O.  So **the `fraction` of each element sums to 1** — an identity that can be checked the moment the file is opened.

Each element appears under **one characterisation only**: Qn for Qn network formers, coordination number
for the other formers and the modifiers, type for ligands.  P has no coordination row here (that is in `network_coordination.csv`).

The `O_b` row is aggregated from three rows of `ligand_type` (Al-Al / Al-P / P-P), and its `sd` is
**re-accumulated frame by frame** rather than summed over the three — the variance of a sum of correlated terms is not the sum of the variances.

The former `network_average.csv` has been removed.  Its two means are still in the `[inputs]` block of every file
(`mean_qn P=0.95  mean_n_bo Al=4.12 P=2.40  mean_cn Al=4.12 P=4.00 Zn=4.10`),
and can also be recomputed exactly from this table: $\sum_n n \cdot f_n = 0.946$.  Note that `mean_qn`
(homonuclear connections) and `mean_n_bo` (bridging-oxygen count) are two different quantities; they are given side by side precisely to make the statement obvious.

### Example: `network_qn.csv`

```csv
file,label,former,qn,count,fraction,sd
43Z43P15A,P-Q0,P,0,465,2.500000e-1,0.000000e0
43Z43P15A,P-Q1,P,1,1035,5.564516e-1,0.000000e0
43Z43P15A,P-Q2,P,2,355,1.908602e-1,0.000000e0
43Z43P15A,P-Q3,P,3,5,2.688172e-3,0.000000e0
```

Al is not among them — in the same run it appears in `network_coordination.csv` on the `cn` 4/5/6 rows.

### Example: `network_ligand_type.csv`

```csv
file,label,type,former_a,former_b,count,fraction,sd
43Z43P15A,P-O_n,O_n,P,,2985,4.543379e-1,1.203302e-3
43Z43P15A,Al-O_b-Al,O_b,Al,Al,10,1.522070e-3,0.000000e0
43Z43P15A,Al-O_b-P,O_b,Al,P,2693,4.098935e-1,1.154167e-3
43Z43P15A,P-O_b-P,O_b,P,P,880,1.339422e-1,0.000000e0
43Z43P15A,O_t,O_t,,,2,3.044140e-4,4.168360e-4
```

`label` reads the whole row as one phrase `<former_a>-<type>-<former_b>`; missing positions get **no
placeholder**, so a free oxygen is just `O_f`, a non-bridging one is `P-O_n` and a tricluster is `O_t`.
`former_a` / `former_b` remain separate columns, so querying Al-O-Al is still
`query("former_a=='Al' and former_b=='Al'")` without slicing strings.

The partner columns of `O_t` are left empty — it is not a pair, and its pairings are given by the
`n_formers=3` rows of the `linkage` table.  Note there is **no `Al-O_n` row** here: the Al in this system
carries no non-bridging oxygen, which is exactly why its number of bridges and its coordination number agree everywhere.

> **The denominator of `fraction` is the atom count of that ligand element**, not of all ligand atoms.  This
> is the statement behind the BO fraction in the literature; "bridging oxygens as a fraction of O+F"
> corresponds to no commonly used quantity.  The two coincide with a single ligand species, and differ only in systems like `--Al-O` + `--Al-F`.

### Example: `network_linkage.csv`

```csv
file,linkage,ligand,elem_a,qn_a,cn_a,elem_b,qn_b,cn_b,n_formers,count,...
43Z43P15A,Al_4-O-P_0,O,Al,0,4,P,0,4,2,844,...
43Z43P15A,Al_4-O-P_1,O,Al,0,4,P,1,4,2,1210,...
43Z43P15A,Al_4-O-P_2,O,Al,0,4,P,2,4,2,213,...
```

The 4 in `Al_4` is the **coordination number** (read `cn_a`), the 1 in `P_1` is **$n$, the P–O–P count of
that P** (read `qn_b`) — this is the literature's own convention (Al[4] against $Q^n$), and the `#` header
of each file says so.  Rows like `Al_4-O-P_0` are common: that P is connected to an Al through this bridge,
but has no P–O–P of its own and is therefore $Q^0$.  A `qn_a=0` on the Al end means there is no Al–O–Al.

The **atom vocabulary** is used here (`P_3`, not `P-Q3`): a linkage describes the connection between two
atoms, while Qn names a structural unit containing several atoms.  This matches the labels of the exported trajectory.

### Analysing `linkage` with pandas

```python
import pandas as pd
d = pd.read_csv("network_linkage.csv", comment="#")

# aggregate by element pair
d.groupby(["elem_a", "elem_b"])["count"].sum()

# which Al coordinations Al-O-Al occurs between — directly comparable to Al[4]/Al[5]/Al[6] in NMR
al = d.query("elem_a == 'Al' and elem_b == 'Al' and n_formers == 2")
al.pivot_table(index="cn_a", columns="cn_b", values="count", aggfunc="sum")

# Qn-Qn connection matrix of P-O-P
d.query("elem_a == 'P' and elem_b == 'P'").pivot_table(
    index="qn_a", columns="qn_b", values="count", aggfunc="sum")

# which Qn of P an Al of a given coordination prefers — compare with 27Al Al[4]/[5]/[6] x 31P Qn
d.query("elem_a == 'Al' and elem_b == 'P'").pivot_table(
    index="cn_a", columns="qn_b", values="count", aggfunc="sum")

# true bridges only, excluding three-coordinate ligands
d.query("n_formers == 2")

# several ligand species: keep O bridges and F bridges apart
d.groupby("ligand")["count"].sum()
```

> **Cross-checks.** The `n_formers=2` count of `Al-O-Al` in `linkage` should equal the `O_b, Al, Al` row of
> the `ligand_type` table; the difference comes from three-coordinate oxygens.
> $\sum(m_\mathrm{Al} \times \text{count})$ (from `qn_partner`) should equal the true-bridge count of
> `Al-O-P` in `linkage`.  When neither holds, check first whether a cutoff was left out.

---

## Labelled trajectory (`--export-traj`)

One file per input: `<input stem>_types[_<suffix>].<ext>`.  The input stem has to be part of the name,
otherwise in batch mode the second input would overwrite the first (the same reason as for the cubes of `ferro map`).

### Labels

| Role | Label | The number is |
|---|---|---|
| Qn network former | `P_0` `P_1` … `Si_4` | **n** (same-element connection count, P–O–P; heteronuclear connections excluded) |
| other network formers | `Al_4` `Al_5` `Al_6` | **coordination number** |
| free ligand | `O_f` | — |
| non-bridging ligand | `O_n` | — |
| bridging ligand | `O_b` | — |
| three-coordinate ligand | `O_t` | — |
| modifier | `Zn` | no suffix |

`ferro net` prints this table once per run with the actual elements filled in, because which element
reports Qn and which reports coordination number depends on `--qn` and on which cutoffs were given.

### Two vocabularies: unit and atom

| Vocabulary | Looks like | Used in | Why |
|---|---|---|---|
| **unit** | `P-Q2`, `Al_4` | `composition` / `qn` / `qn_partner` | these three tables count "how many Q2 units there are", and Qn is by origin the notation for a structural unit |
| **atom** | `P_2`, `Al_4` | `linkage`, the exported trajectory | a linkage describes a connection **between atoms**, while one Qn unit contains several atoms |

For non-Qn network formers the two vocabularies happen to coincide (both `Al_4`), because Al has no unit notation available.

> **Downstream type selection uses the atom vocabulary.** In the exported trajectory that P is called
> `P_2`, not `P-Q2` — a trajectory label has to be splittable back into an element by the dump reader at
> the first underscore, and the `Q` of `Q2` is not an element.  So it is `traj gr -x P_2`, not `-x P-Q2`.

The element prefix (`P-Q2` rather than `Q2`) is there because in systems with two Qn network formers
(B+P in borophosphates, Si+Al in aluminosilicates) `Q0`…`Q4` would each appear twice, and **no neighbouring column can resolve a name collision across elements**.

All of them follow the `<element>_<suffix>` convention and are split at the **first underscore**.  Modifiers
carry no role suffix: the old scheme binned them by non-bridging-ligand count 0/1/2/≥3, while the actual
coordination number of a modifier is 3–6, and measured on the reference trajectory 97% fell into the fallback bin — no resolving power at all.

> **The same-looking number means different things per element.** This is the literature's own reading
> ($Q^2$ and Al[4] are never confused), at the price that in a system where Al carries no non-bridging
> oxygen the two statements give the same number and a mistake would be invisible — which is why the convention is written into the type at classification time rather than looked up at rendering time: change the list with `--qn` and the labels change with it.

### How the two formats differ

| | `lammpstrj` (default) | `extxyz` |
|---|---|---|
| where the label lives | folded into the `element` column | its own `label:S:1` column |
| `species` / `element` | the label (`P_2`) | the bare element (`P`) |
| reading back | the reader splits it at the first underscore into element + label | read straight from separate columns, no guessing |
| downstream compatibility | no new column, existing tools parse it as before | one new column |

The dump folds because it has only one name column to put a second name in; the columns of extxyz are
self-describing and need no such compromise.  **Only `ferro net --export-traj` folds** — `ferro convert`
writes clean element symbols whatever the labels are.  The fold also has a guard: a label not of the form
`<element>_…` is not folded, and those occurrences are counted and warned about (a CIF `O1` or a CP2K `Fe1` folded in and read back would be a non-existent element `O1`).

### Downstream type selection

```bash
# select by element — works with any number of frames
ferro traj gr -i run_types.lammpstrj -a P -b O

# select by label — requires a single frame
ferro traj gr -i run_types.lammpstrj -x P_3 -y O_b --last-n 1

# bridging oxygens around five-coordinate Al — the entry point for locating Al-O-Al
ferro traj gr -i run_types.lammpstrj -x Al_5 -y O_b --last-n 1
```

> **Selecting by label only holds for a single frame.** `g(r)` requires the per-type particle count to be
> conserved, but labels are dynamic — a P changes its own Qn as the trajectory advances (measured on the
> reference trajectory, `P_3` is 149 / 152 / 150 / 150 / 150 frame by frame).  Selecting by label on a
> multi-frame labelled trajectory is rejected with a "per-type atom counts change" error; that is a guard, not a defect.  Select by element for multiple frames.

The labelled trajectory and the statistics tables are independent outputs and can verify each other:

- the first-shell CN of `traj gr -x P_3 -y O_b --last-n 1` should be **exactly 3** — the classification layer
  says "this P has 3 bridging oxygens", and the independent CN integration path should count 3 back.
- `traj gr -x Al_5 -y O_b --last-n 1` measures 4.933 rather than 5.000 on the reference trajectory: the
  missing 0.067 is one three-coordinate oxygen, labelled `O_t` and not `O_b`.  It adds up, which confirms that the number in the label really is the coordination number.

---

## Typical parameters for reference

| System | Network former | Modifier | Reference cutoff |
|---|---|---|---|
| $P_2O_5$ glass | P | — | P-O: 2.3–2.4 Å |
| $SiO_2$ glass | Si | — | Si-O: 1.8 Å |
| $Al_2O_3$ | Al | — | Al-O: 2.1–2.4 Å |
| $GeO_2$ glass | Ge | — | Ge-O: 2.0 Å |
| ZnO–$P_2O_5$ glass | P | Zn | P-O: 2.4 Å, Zn-O: 2.6 Å |
| ZnO–$Al_2O_3$–$P_2O_5$ | P, Al | Zn | P-O: 2.4, Al-O: 2.4, Zn-O: 2.6 Å |

Cutoffs should be chosen from the position of the first minimum in g(r):

```bash
ferro traj gr -i traj.lammpstrj -a P -b O --r-max 5 --plot
```
