# ferro-py Configuration

`ferro-py job` and `ferro-py conv` need to know where the pseudopotentials,
CP2K's data files and your own templates are. That is the one global
configuration file:

```
~/.config/ferro/config.toml          # or the path in $FERRO_CONFIG
```

```bash
ferro-py config init                 # write a commented template there
ferro-py config init -o cfg.toml     # somewhere else (then export FERRO_CONFIG=cfg.toml)
```

`config init` never overwrites an existing file without `--force`; `--mkdir`
creates a missing directory for `-o` without asking.

## Every entry is optional

A missing file is an empty configuration. An entry is needed only when a job
uses it, and the error then names the entry to add — a VASP job without
`[vasp] potcar_dir` stops with

```
[vasp] potcar_dir is not set in /home/me/.config/ferro/config.toml (create the file with `ferro-py config init`)
```

A misspelt key (`potcar_dri`) or an unknown section is an error, not ignored:
otherwise the file looks configured while the default is silently in use.

Paths accept `~` and `$VARS`; a relative path is relative to the configuration
file itself, so a project-local configuration can say `potcar_dir = "pp"`.

## Entries

| Section | Key | Meaning |
|---|---|---|
| `[vasp]` | `potcar_dir` | One sub-directory per POTCAR: `<potcar_dir>/<name>/POTCAR` |
| `[vasp.potcar]` | `<Element> = "<name>"` | The sub-directory for one element, overriding the recommended one |
| `[qe]` | `pseudo_dir` | The directory holding the UPF files; written into pw.x's `pseudo_dir` |
| `[qe.pseudos]` | `<Element> = "<file>"` | The UPF file for an element that has several in `pseudo_dir` |
| `[cp2k]` | `data_dir` | CP2K's `data/` (`BASIS_MOLOPT`, `GTH_POTENTIALS`, `BASIS_pob`, `ALL_POTENTIALS`, …) |
| `[cp2k]` | `basis_family` | Pseudopotential runs: basis family  [`DZVP-MOLOPT-SR-GTH`] |
| `[cp2k]` | `potential_family` | Pseudopotential runs: GTH family  [`GTH-PBE`] |
| `[cp2k]` | `ae_basis_family` | All-electron runs: basis family  [`pob-TZVP-rev2`] |
| `[cp2k.kinds]` | `<Element> = { basis = "…", potential = "…" }` | Basis set and / or potential of one element, by name |
| `[templates]` | `dir` | Your own templates (below) |

How each program uses its entries — which POTCAR a VASP job takes, how the UPF
files are indexed, how CP2K pairs a basis set with a potential — is on the
program's page: `ferro-py doc job vasp`, `job qe`, `job cp2k`.

## Your own templates

A template is a native input file of the program (see `ferro-py doc job`). Put
your edited ones under `[templates] dir`, laid out as

```
<dir>/vasp/{sp,opt,aimd}/INCAR
<dir>/qe/{sp,opt,aimd}/pw.in
<dir>/cp2k/{sp,opt,aimd}/cp2k.inp
```

and a file there replaces the built-in template of the same program and task for
every job; the ones you do not provide stay built-in. `--template FILE` on the
command line beats both, for that one run.

## Example

```toml
[vasp]
potcar_dir = "~/pp/vasp/potpaw_PBE.64"

[vasp.potcar]
Zn = "Zn"            # instead of the recommended one

[qe]
pseudo_dir = "~/pp/qe/pslibrary"

[qe.pseudos]
O = "O.pbe-n-kjpaw_psl.1.0.0.UPF"

[cp2k]
data_dir = "~/soft/cp2k-2026.2/data"
potential_family = "GTH-PBE"

[cp2k.kinds]
Zn = { basis = "DZVP-MOLOPT-SR-GTH-q12", potential = "GTH-PBE-q12" }

[templates]
dir = "~/.config/ferro/templates"
```
