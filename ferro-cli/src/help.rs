
pub fn print_fe_job_overview() {
    println!(
        r#"ferro job — Generate QC software input files

Usage:
  ferro job -s <SOFTWARE> -i <FILE> [OPTIONS]
  ferro job -s <SOFTWARE>              show software-specific parameters

Supported software:
  gaussian   Gaussian 16/09 input file (.gjf)
  cp2k       CP2K input file (.inp)  — DFT/MD/GeoOpt/CellOpt
  qe         Quantum ESPRESSO pw.x input (.in)  — scf/relax/md/bands

Common options:
  -i, --input  PATH   Input structure file (xyz, cif, pdb, POSCAR, …)
  -o, --output PATH   Output file (default: job.gjf / job.inp)
      --metal-units   LAMMPS metal units for dump files

Full documentation:  ferro doc job"#
    );
}

pub fn print_job_help(software: &str) {
    match software.to_lowercase().as_str() {
        "gaussian"                     => print_job_gaussian(),
        "cp2k"                         => print_job_cp2k(),
        "qe" | "espresso" | "pwscf"    => print_job_qe(),
        other => println!("Unknown software: {other}  (supported: gaussian | cp2k | qe)"),
    }
}

fn print_job_qe() {
    println!(
        r#"ferro job -s qe — Quantum ESPRESSO pw.x input file

  ibrav = 0; cell taken from the structure (CELL_PARAMETERS angstrom).
  Pseudopotentials are referenced as <Element>.UPF in --pseudo-dir.

Task and electronic structure:
  --qe-task STR       scf | nscf | bands | relax | vc-relax | md | vc-md  [scf]
  --qe-functional STR pbe pbesol revpbe blyp scan r2scan pbe0 hse06       [pbe]
  --ecutwfc F         Plane-wave cutoff [Ry]                               [50]
  --smearing STR      none | gaussian | mp | mv | fd                      [none]
                      (mp / mv are the ones to use for metals)
  --kpoints K1 K2 K3  Monkhorst-Pack mesh (omit for Gamma)
  --pseudo-dir PATH   Pseudopotential directory                      [./pseudo]

Charge / spin (shared by all three targets):
  --charge INT        Override total charge
  --multiplicity INT  Override 2S+1 (-> nspin=2, tot_magnetization)
  --auto-spin         Guess it from the structure; ON by default for qe

MD (--qe-task md|vc-md):
  --md-steps INT      Number of MD steps                               [10000]
  --temperature F     Target temperature [K]                          [298.15]

Examples:
  ferro job -s qe -i crystal.cif
  ferro job -s qe -i metal.cif --smearing mp --kpoints 8 8 8
  ferro job -s qe -i slab.xyz --qe-task relax --qe-functional scan
  ferro job -s qe -i Fe2O3.cif --auto-spin --kpoints 4 4 4 -o pw.in

Full documentation:  ferro doc job"#
    );
}

fn print_job_gaussian() {
    println!(
        r#"ferro job -s gaussian — Gaussian 16/09 input file

Parameters:
  -m, --method  STR   DFT functional           default: B3LYP
  -b, --basis   STR   Basis set                default: 6-31G*
  -o PATH             Output file              default: job.gjf

Charge / spin (shared):
  --charge INT        Override total system charge
  --multiplicity INT  Override spin multiplicity 2S+1 (highest priority)
  --auto-spin         Guess multiplicity from structure:
                        magmom 求和 → 氧化态+Hund → 电子数奇偶下限

Example:
  ferro job -s gaussian -i mol.xyz
  ferro job -s gaussian -i mol.xyz -m PBE0 -b def2-TZVP -o sp.gjf
  ferro job -s gaussian -i FeCl3.xyz --auto-spin            # 推断高自旋多重度
  ferro job -s gaussian -i radical.xyz --charge 0 --multiplicity 2

Full documentation:  ferro doc job"#
    );
}

fn print_job_cp2k() {
    println!(
        r#"ferro job -s cp2k — CP2K input file (GPW/DFT, periodic systems)

Task and electronic structure:
  --task STR          energy | force | geo-opt | cell-opt | md | freq  [energy]
  --functional STR    pbe blyp revpbe pbesol          (GGA)             [pbe]
                      pbe0 b3lyp hse06                (hybrid, auto &HF block)
                      scan r2scan                     (meta-GGA via LIBXC)
  --cp2k-basis STR    dzvp-molopt-sr tzvp-molopt tzv2p-molopt   [dzvp-molopt-sr]
                      dzvp-gth tzvp-gth               (older GTH style)
                      pob-dzvp pob-tzvp               (all-electron, periodic)
  --dispersion STR    none | d3 | d3bj                                  [none]
  --scf STR           diag (metals, large systems) | ot (band-gap systems) [diag]
  --cutoff INT        Plane-wave cutoff [Ry]                              [400]
  --rel-cutoff INT    Relative cutoff [Ry]                                 [50]
  --smear             Fermi-Dirac smearing (300 K)
  --pbc STR           xyz | z | none                              (auto from cell)
  --kpoints K1 K2 K3  Monkhorst-Pack mesh

Charge / spin (shared by all three targets):
  --charge INT        Override total system charge
  --multiplicity INT  Override 2S+1 (highest priority; disables auto-spin)
  --auto-spin         Guess it from the structure (see `ferro doc spin`)

Output:
  --atom-charge STR   none | mulliken | hirshfeld | hirshfeld-i         [none]
  --cube STR          none | density | elf | hartree                    [none]
  --molden            Export a Molden wavefunction file
  --project STR       CP2K project name                                [ferro]

MD (--task md):
  --md-steps INT      Number of MD steps                               [10000]
  --md-timestep F     Timestep [fs]                                       [1.0]
  --temperature F     Temperature [K]                                  [298.15]
  --thermostat STR    csvr | nose | langevin | none (NVE)                [csvr]
  --traj-freq INT     Write the trajectory every N steps                  [100]
  --barostat          NPT with a flexible cell

  The exact per-element basis and pseudopotential names are resolved from a
  2829-entry database; --cp2k-basis only picks the family.

Examples:
  ferro job -s cp2k -i glass.xyz
  ferro job -s cp2k -i glass.xyz --task geo-opt --dispersion d3bj -o opt.inp
  ferro job -s cp2k -i glass.xyz --task md --temperature 1500 --md-steps 50000
  ferro job -s cp2k -i crystal.cif --functional pbe0 --scf ot --cp2k-basis pob-tzvp

Full documentation:  ferro doc job"#
    );
}

/// `ferro convert` with no `-i`: the read/write format matrix.
pub fn print_convert() {
    println!(
        r#"ferro convert — Structure / trajectory format conversion

  Reads one file, writes another. Both formats come from the file NAMES;
  there is no --from / --to flag.

Supported formats:
{}

Parameters:
  -i, --input  FILE       Input file  (format from its name)
  -o, --output FILE       Output file — a full PATH here, unlike the analysis
                          commands where -o is a suffix
      --start  N          First frame to take      (0-based, inclusive) [0]
      --end    N          Last frame to take       (0-based, INCLUSIVE) [last]
      --stride N          Take every Nth frame within [start, end]      [1]
      --number N          Take this many frames, spread evenly, both ends kept
      --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces
                          eV/Å); default is real units
  -h, --help              Short parameter table (this page adds the formats)

Frame selection:
  --end is INCLUSIVE and 0-based, matching the numbers `ferro info` prints.
  --stride (a spacing) and --number (a total) cannot be combined.

How many files come out — decided by the target format, not by a flag:
  one file          if the format holds a trajectory (Frames column above)
  one file PER FRAME if it holds a single structure (POSCAR, data, QE)
                    -> POSCAR_0000, conf_0002.vasp; the number is the index in
                       the ORIGINAL trajectory, zero-padded to 4 digits

Examples:
  ferro convert -i input.cif -o POSCAR
  ferro convert -i traj.dump -o sub.extxyz --start 100 --end 199
  ferro convert -i traj.dump -o conf.vasp --number 20

Full documentation:  ferro doc convert"#,
        crate::io_dispatch::supported_formats()
    );
}

/// `ferro info` with no `-i`: what the summary reports.
pub fn print_info() {
    println!(
        r#"ferro info — Structure / trajectory summary

  Frame count, composition, cell parameters, volume and mass density.
  Reads every format `ferro convert` reads (run `ferro convert` for the list).

Parameters:
  -i, --input  FILE       Input file (format from its name)
      --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces
                          eV/Å); default is real units
  -h, --help              Short parameter table

Reported for the FIRST and the LAST frame only, not every frame:
  atoms + per-element composition · cell a b c / α β γ · volume · density
  (g/cm³) · per-axis PBC flags · whether energy / forces / velocities are there

  Density is omitted when there is no cell — undefined without a volume, and a
  placeholder would read like a measurement. For a mean ± σ over ALL frames,
  read the `# volume = <mean> +/- <std>` header any `ferro traj` product carries.

  A symbol missing from the element table falls back to 1 amu, which drags the
  density DOWN with no other visible sign — so a warning naming how many atoms
  fell back follows the density line. Treat the number as invalid until it goes.

Examples:
  ferro info -i input.xyz
  ferro info -i traj.lammpstrj

Full documentation:  ferro doc info"#
    );
}

/// `ferro bader` with no `-i`: methods, outputs, and the file-name collision.
pub fn print_bader() {
    println!(
        r#"ferro bader — Bader charge partitioning from a DFT charge density

  Partitions the charge density into atomic basins along its zero-flux surfaces
  and reports the charge, volume and surface distance of each.

Input (format from the file name):  .cube = Gaussian / QE pp.x   anything else = VASP CHGCAR

Parameters:
  -i, --input  FILE       Charge density file
  -m, --method NAME       ongrid | neargrid | offgrid | weight     [neargrid]
  -r, --refine INT        Edge refinement: -1 auto, -2 single pass, N passes [-1]
  -v, --vacval FLOAT      Vacuum density threshold [e/Å³]             [1e-3]
  -h, --help              Short parameter table

Choosing a method:
  neargrid   Default; accurate on ordinary cells
  ongrid     Cheapest, but staircased basin surfaces bias the charges
  offgrid    Interpolated gradients — slower, no grid bias
  weight     Yu-Trinkle; use it for strongly skewed (non-orthogonal) cells

Output — three Henkelman-format .dat files named after the INPUT stem
  <stem>_ACF.dat  <stem>_BCF.dat  <stem>_AVF.dat
  Kept in the Henkelman layout (not csv) because external tools parse them.

  CAUTION: written to the CURRENT DIRECTORY; there is no --outdir yet. VASP
  names every charge density CHGCAR, so two systems run from one working
  directory both write CHGCAR_ACF.dat and the second silently overwrites the
  first. Until --outdir lands, cd into each system's directory.

Examples:
  ferro bader -i CHGCAR
  ferro bader -i CHGCAR --method weight
  ferro bader -i CHGCAR --method neargrid --refine 3 --vacval 1e-4

Full documentation:  ferro doc bader"#
    );
}

/// Top-level `ferro` help: the command families, grouped by what they produce.
pub fn print_overview() {
    println!(
        r#"ferro — Computational Chemistry Toolkit  v{}

Usage:
  ferro <GROUP> <COMMAND> [OPTIONS]
  ferro <GROUP>                    list that group's commands
  ferro <GROUP> <COMMAND>          show that command's parameters

Trajectory analysis      one stacked csv per run, `file` as a column; --plot for a look
  traj gr        Radial distribution g(r) + coordination number CN(r)
  traj sq        Structure factor S(q)
  traj msd       Mean square displacement + self-diffusion D
  traj angle     Bond angle distribution P(theta)
  traj vacf      Velocity autocorrelation + Green-Kubo diffusion
  traj rotcorr   Rotational correlation C2(t)
  traj vanhove   Van Hove self-correlation Gs(r,tau)

Spatial maps             one .cube grid file per input — no summary table, no plot
  map density | velocity | force | radius | sdf | chg-sdf

Topology & charges
  net            Structural composition, Qn speciation, ligand types,
                 coordination numbers, bridge connectivity;
                 --export-traj also writes the classified trajectory
  bader          Bader charge partitioning (ACF/BCF/AVF, Henkelman format)

Structure I/O
  convert        Format conversion — run it bare for the read/write matrix
  info           Atoms, cell, volume, density (g/cm³) of a structure or trajectory
  job            Quantum-chemistry input files (gaussian | cp2k | qe)

Machine-learning datasets
  dataset collect   AIMD output -> DeePMD system directories (set.*/*.npy)
  dataset filter    Drop low-quality frames (force / stress thresholds)
  dataset merge     Combine same-composition datasets, shuffle, resize sets

Manual
  doc            The user manual, in the binary. `ferro doc` lists the topics;
                 `ferro doc dataset filter` reads one.

Batch input:
  -i takes several files and expands glob patterns itself — quote them:
    ferro traj gr -i 'runs/*/prod.dump' -a P -b O -o scan
  Each input is analysed on its own; results stack into ONE csv with a `file`
  column. A failed input is skipped, its reason printed, exit code set to 1.

Output naming:
  <outdir>/<command>[_<table>][_<label>]_<suffix>.csv
  --outdir DIR   where every product goes (created if missing; default: cwd)
  -o SUFFIX      batch tag, chosen by you
  <label>        what was analysed: `traj gr -a P -b O` -> gr_P-O.csv, no
                 selection -> gr_all.csv. Label before suffix, so
                 `ls gr_P-O_*` lists one pair across every batch.

  `dataset` is the exception: its products are DIRECTORIES, so -o is an output
  root rather than a suffix.

Help:
  A command typed without -i prints its own page; `-h` gives the short parameter
  table; `ferro doc <topic>` gives the full manual page.

Full documentation:  ferro doc cli-reference"#,
        env!("CARGO_PKG_VERSION")
    );
}

/// `ferro traj` with no subcommand.
pub fn print_traj_overview() {
    println!(
        r#"ferro traj — Trajectory analysis

Usage:
  ferro traj <COMMAND> -i <FILE> [FILE ...] [OPTIONS]
  ferro traj <COMMAND>             show that command's parameters

Commands:
  gr        Radial distribution function g(r) and coordination number CN(r)
  sq        Structure factor S(q) via Fourier transform of g(r)
  msd       Mean square displacement MSD(t), time-shift averaged
  angle     Bond angle distribution P(theta) for A-B-C triplets
  vacf      Velocity autocorrelation function + Green-Kubo diffusion
  rotcorr   Rotational correlation C2(t) for molecular bond vectors
  vanhove   Van Hove self-correlation Gs(r, tau)

Common options:
  -i, --input  FILE...  Input trajectory file(s); glob patterns allowed (quote them)
  -o, --output SUFFIX   Batch tag -> <command>[_<table>][_<label>]_<suffix>.csv
      --outdir DIR      Write every product here (created if missing; default: .)
      --last-n N        Use only the last N frames of the trajectory
      --ncore  N        Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units (velocities in A/ps)

Selecting types (gr / angle only):
  -a -b -c              by chemical element   (e.g. -a P -b O)
  -x -y -z              by site label         (e.g. -x P_2 -y O_b)
  The two groups are mutually exclusive. For gr the first is the centre and the second
  the neighbour, and the order matters for CN. sq has no type selection: its partials
  sum back to the weighted totals, so filtering to one pair hides that closure."#
    );
}

/// `ferro map` with no subcommand.
pub fn print_map_overview() {
    println!(
        r#"ferro map — Spatial distribution maps

Usage:
  ferro map <COMMAND> -i <FILE> [FILE ...] [OPTIONS]
  ferro map <COMMAND>              show that command's parameters

Commands:
  density   Time-averaged number density [atoms/A^3] per voxel
  velocity  Time-averaged speed |v| per voxel  (needs frame velocities)
  force     Time-averaged force magnitude |f| per voxel  (needs frame forces)
  radius    Hard-sphere spatial occupancy map
  sdf       Cluster spatial distribution function (Qn-type, Kabsch alignment)
  chg-sdf   Charge-density cluster SDF from QE pp.x cube files (--cubes, not -i)

Common options:
  -i, --input  FILE...  Input trajectory file(s); glob patterns allowed (quote them)
  -o, --output STEM     Output file stem (default depends on the command)
      --last-n N        Use only the last N frames
      --ncore  N        Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units (velocities in A/ps)

Multiple inputs:
  Unlike ferro traj, the product is one 3-D grid file per input — there is nothing to
  stack. File names therefore carry the input stem (density_<stem>.cube) so several
  inputs cannot overwrite each other. No summary table, no plot."#
    );
}

// ─── ferro traj: gr / sq / msd / angle ──────────────────────────────────────

pub fn print_gr() {
    println!(
        r#"ferro traj gr — Radial Distribution Function and Coordination Number

  g(r) and CN(r) for every ordered pair of types, into one file.
  Requires a periodic cell.

  g(r) is symmetric (A-B = B-A); CN(r) is DIRECTED — A-B is the average number
  of B around each A, so the two generally differ.

Selecting a pair (centre first, neighbour second):
  -a ELEM  -b ELEM        by element:    -a P -b O  -> O around each P
  -x LABEL -y LABEL       by site label:  -x P_2 -y O_b
  (mutually exclusive; omit both for every pair)

  NOTE: -x/-y hold over a SINGLE frame only — g(r) needs a fixed particle count
  per type and labels shift as the run evolves. Use --last-n 1, or select by
  element.

Parameters:
  --r-min  FLOAT          Min cutoff radius [Å]                     [0.001]
  --r-max  FLOAT          Max cutoff radius [Å]                    [10.005]
                          (clamped to half the smallest interplanar spacing)
  --dr     FLOAT          Histogram bin width [Å]                   [0.002]
  --last-n INT            Use only the last N frames
  --ncore  INT            Parallel threads                    [all cores]
  -o SUFFIX               Batch tag  -> gr_<pair>_<suffix>.csv
  --outdir DIR            Write products here (created if missing)
  --metal-units           LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)
  --plot                  PNG next to the data file (needs a pair)

Output — long format, one row per (file, pair, r): file pair r g_r cn_r
  File name: gr_<pair>[_<suffix>].csv; no selection -> gr_all.csv

Examples:
  ferro traj gr -i traj.dump -a P -b O
  ferro traj gr -i 'runs/*/prod.dump' -a P -b O -o scan
  ferro traj gr -i traj.dump -x Al_5 -y O_b --last-n 1

Full documentation:  ferro doc traj gr"#
    );
}

pub fn print_sq() {
    println!(
        r#"ferro traj sq — Structure Factor S(q)

  S(q) by Fourier transform of g(r) (Faber-Ziman), weighted by XRD
  (Waasmaier-Kirfel) form factors or neutron scattering lengths.

  No type selection: EVERY pair is written, always. The primary product is the
  pair of totals; the partials are a decomposition that sums back to them, and
  keeping one pair would hide exactly that closure. Filter columns in pandas.

Parameters:
  --q-min      FLOAT      Min q [Å⁻¹]                                 [0.1]
  --q-max      FLOAT      Max q [Å⁻¹]                                [25.0]
  --dq         FLOAT      q bin width [Å⁻¹]                          [0.02]
  --weighting  ENUM       none | xrd | neutron | both                [both]
  --r-min      FLOAT      g(r) lower cutoff [Å]                     [0.001]
  --r-max      FLOAT      g(r) cutoff [Å]                          [10.005]
  --dr         FLOAT      g(r) bin width [Å]                        [0.002]
  --last-n     INT        Use only the last N frames
  --ncore      INT        Parallel threads (used in the g(r) step)
  -o SUFFIX               Batch tag -> sq_<suffix>.csv
  --outdir DIR            Write products here (created if missing)
  --metal-units           LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)
  --plot                  PNG next to the data file (weighted totals only)

Output — WIDE format, one row per (file, q):
  file  q  total_xrd  total_neutron  then three columns per pair
  (_sq unweighted, _xrd and _neutron carrying w_ij(q)*S_ij(q))

  Wide, not long like gr: the totals are one value per q. Inputs with different
  element sets contribute different pair columns; gaps stay empty (NaN).

Examples:
  ferro traj sq -i traj.dump
  ferro traj sq -i traj.dump --weighting xrd --q-max 20.0 -o xrd
  ferro traj sq -i 'runs/*/prod.dump' -o scan

Full documentation:  ferro doc traj sq"#
    );
}

pub fn print_msd() {
    println!(
        r#"ferro traj msd — Mean Square Displacement
  Computes MSD(t) = <|r(t₀+t) − r(t₀)|²> averaged over time origins.
  Outputs total MSD and per-axis (a/b/c) components.

Parameters:
  --dt        FLOAT      Timestep between frames [fs]   default: 1.0
  --shift     INT        Time-origin stride             default: 1
  --elements  Fe,O,...   Track only these elements      default: all
  --fit-range FMIN,FMAX  Linear-fit window as fractions of the MSD
                         curve; reports self-diffusion D = slope/6
                         (Einstein, 3-D) and R²
  --last-n    INT        Use only the last N frames
  --ncore     INT        Parallel threads
  --plot                 Generate PNG and open in viewer
  -o SUFFIX              Batch tag -> msd_<elements>_<suffix>.csv
  --outdir DIR           Write products here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

File name — msd_<elements>[_<suffix>].csv, elements sorted:
  --elements P,O -> msd_O-P.csv     (so does --elements O,P: it is a set, and the
                                     same data must not land under two names)
  no --elements  -> msd_all.csv

Example:
  ferro traj msd -i traj.xyz --dt 2.0
  ferro traj msd -i traj.dump --elements Li --dt 1.0 --last-n 2000
  ferro traj msd -i traj.dump --dt 1.0 --fit-range 0.3,0.8 --plot

Full documentation:  ferro doc traj msd"#
    );
}

pub fn print_angle() {
    println!(
        r#"ferro traj angle — Bond Angle Distribution

  P(θ) for all A-B-C triplets within cutoff distances. B is the central atom.

Selecting a triplet (all three required):
  -a ELEM  -b ELEM  -c ELEM      by element,    B is the centre
  -x LABEL -y LABEL -z LABEL     by site label, Y is the centre
  (mutually exclusive)

Parameters:
  --r-cut-ab  FLOAT       End-A-to-centre-B cutoff [Å]                 [2.3]
  --r-cut-bc  FLOAT       End-C-to-centre-B cutoff [Å]                 [2.3]
                          A is what -a/-x names, C what -c/-z names. Without a
                          named triplet both fall back to canonical (Z, symbol)
                          order; equal end types take min(ab, bc).
  --angle-min FLOAT       Histogram lower edge [°]                     [0.0]
  --angle-max FLOAT       Histogram upper edge [°]                   [180.0]
                          Angles outside the window are DISCARDED, not hidden.
  --d-angle   FLOAT       Histogram bin width [°]                      [0.1]
  --last-n    INT         Use only the last N frames
  --ncore     INT         Parallel threads                       [all cores]
  -o SUFFIX               Batch tag -> angle_<triplet>_<suffix>.csv
  --outdir DIR            Write products here (created if missing)
  --metal-units           LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)
  --plot                  PNG next to the data file

Output — long format: file triplet theta count p
  File name: angle_<triplet>[_<suffix>].csv, triplet as you wrote it
  (-a O -b P -c O -> angle_O-P-O.csv); no triplet -> angle_all.csv

  Each geometric angle is counted ONCE: a PO4 tetrahedron gives 6 O-P-O angles,
  not 12.

Examples:
  ferro traj angle -i traj.dump -a O -b P -c O --r-cut-ab 2.0 --r-cut-bc 2.0
  ferro traj angle -i traj.dump -x O_b -y P_2 -z O_n --last-n 1
  ferro traj angle -i traj.dump -a O -b P -c O --angle-min 90 --angle-max 130

Full documentation:  ferro doc traj angle"#
    );
}

// ─── ferro traj: vacf / rotcorr / vanhove ───────────────────────────────────

pub fn print_vacf() {
    println!(
        r#"ferro traj vacf — Velocity Autocorrelation Function
  Computes C_v(t) = <v(t₀)·v(t₀+t)> / <v²(t₀)>, averaged over origins.
  Also outputs running integral (Green-Kubo diffusion coefficient).
  Requires frame.velocities in the input file.

Parameters:
  --dt       FLOAT      Timestep [fs]                 default: 1.0
  --shift    INT        Time-origin stride             default: 1
  --elements Fe,O,...   Include only these elements    default: all
  --last-n   INT        Use only the last N frames
  --tau      INT        Lag time in frames             default: half the run
  --ncore    INT        Parallel threads               default: all cores
  -o SUFFIX             Batch tag -> vacf_<elements>_<suffix>.csv
  --outdir DIR          Write products here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

File name — vacf_<elements>[_<suffix>].csv, elements sorted; vacf_all.csv without
  --elements.

Example:
  ferro traj vacf -i traj.dump --dt 2.0
  ferro traj vacf -i traj.dump --elements O --last-n 1000

Full documentation:  ferro doc traj vacf"#
    );
}

pub fn print_rotcorr() {
    println!(
        r#"ferro traj rotcorr — Rotational Correlation Function
  Computes C₂(t) = <P₂(û(t₀)·û(t₀+t))> for molecular bond vectors.
  --center and --neighbor are required to define the bond direction.

Parameters:
  --center    ELEM    Central atom element (required)   e.g. O
  --neighbor  ELEM    Neighbor atom element (required)  e.g. H
  --r-cut     FLOAT   Bond search cutoff [Å]            default: 1.2
  --dt        FLOAT   Timestep [fs]                     default: 1.0
  --shift     INT     Time-origin stride                default: 1
  --last-n    INT     Use only the last N frames
  --tau       INT     Lag time in frames                default: half the run
  --ncore     INT     Parallel threads                  default: all cores
  -o SUFFIX           Batch tag -> rotcorr_<centre>-<neighbour>_<suffix>.csv
  --outdir DIR        Write products here (created if missing)
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

File name — rotcorr_<centre>-<neighbour>[_<suffix>].csv; both are required, so this
  one never falls back to "all". --center O --neighbor H -> rotcorr_O-H.csv

Example:
  ferro traj rotcorr -i traj.xyz --center O --neighbor H
  ferro traj rotcorr -i traj.dump --center O --neighbor H --dt 2.0

Full documentation:  ferro doc traj rotcorr"#
    );
}

pub fn print_vanhove() {
    println!(
        r#"ferro traj vanhove — Van Hove Self-Correlation Function
  Computes Gs(r, τ) = probability distribution of atomic displacements
  over a fixed time lag τ.

Parameters:
  --tau      INT        Lag time in frames              default: half trajectory
  --dt       FLOAT      Timestep [fs]                  default: 1.0
  --shift    INT        Time-origin stride              default: 1
  --r-max    FLOAT      Max displacement [Å]           default: 10.0
  --dr       FLOAT      Bin width [Å]                  default: 0.01
  --elements Fe,O,...   Track only these elements       default: all
  --last-n   INT        Use only the last N frames
  --ncore    INT        Parallel threads                default: all cores
  -o SUFFIX             Batch tag -> vanhove_<elements>_<suffix>.csv
  --outdir DIR          Write products here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

File name — vanhove_<elements>[_<suffix>].csv, elements sorted; vanhove_all.csv
  without --elements.

Example:
  ferro traj vanhove -i traj.xyz --tau 100
  ferro traj vanhove -i traj.dump --elements Li --tau 500 --dt 2.0

Full documentation:  ferro doc traj vanhove"#
    );
}

// ─── ferro map ──────────────────────────────────────────────────────────────

pub fn print_cube_density() {
    println!(
        r#"ferro map density — Spatial Number Density
  Divides the simulation box into nx×ny×nz voxels and computes
  the time-averaged atom number density [atoms/Å³] per voxel.
  Output is a Gaussian cube file (readable by VESTA / VMD).

Parameters:
  --nx INT            Grid points along a axis    default: 50
  --ny INT            Grid points along b axis    default: 50
  --nz INT            Grid points along c axis    default: 50
  --elements Fe,O     Count only these elements   default: all
  --last-n   INT      Use only the last N frames
  --ncore    INT      Parallel threads
  -o STEM             Output name stem            default: density.cube
  --outdir DIR        Write the cubes here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Example:
  ferro map density -i traj.dump
  ferro map density -i traj.dump --nx 100 --ny 100 --nz 100 --elements Li

Full documentation:  ferro doc map density"#
    );
}

pub fn print_cube_velocity() {
    println!(
        r#"ferro map velocity — Spatial Velocity Distribution
  Computes the time-averaged speed |v| per voxel [Å/fs].
  Requires frame.velocities in the input file.

Parameters:
  --nx INT            Grid points along a axis    default: 50
  --ny INT            Grid points along b axis    default: 50
  --nz INT            Grid points along c axis    default: 50
  --elements Fe,O     Include only these elements default: all
  --last-n   INT      Use only the last N frames
  --ncore    INT      Parallel threads
  -o STEM             Output name stem            default: velocity.cube
  --outdir DIR        Write the cubes here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Example:
  ferro map velocity -i traj.dump --nx 80 --ny 80 --nz 80

Full documentation:  ferro doc map velocity"#
    );
}

pub fn print_cube_force() {
    println!(
        r#"ferro map force — Spatial Force Distribution
  Computes the time-averaged force magnitude |f| per voxel [eV/Å].
  Requires frame.forces in the input file.

Parameters:
  --nx INT            Grid points along a axis    default: 50
  --ny INT            Grid points along b axis    default: 50
  --nz INT            Grid points along c axis    default: 50
  --elements Fe,O     Include only these elements default: all
  --last-n   INT      Use only the last N frames
  --ncore    INT      Parallel threads
  -o STEM             Output name stem            default: force.cube
  --outdir DIR        Write the cubes here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Example:
  ferro map force -i traj.dump --elements O

Full documentation:  ferro doc map force"#
    );
}

pub fn print_cube_radius() {
    println!(
        r#"ferro map radius — Hard-Sphere Spatial Occupancy Map
  For each voxel, counts how many (frame, atom) pairs have the selected
  atom within --radius Å of the voxel centre.  Applies the minimum-image
  convention for periodic cells.  Output is a Gaussian cube file.

  Unlike -m density (Gaussian broadening / bin-count), this mode uses a
  hard binary criterion: voxel is marked if any atom overlaps it.

Parameters:
  --nx      INT       Grid points along a axis    default: 50
  --ny      INT       Grid points along b axis    default: 50
  --nz      INT       Grid points along c axis    default: 50
  --radius  FLOAT     Hard-sphere cutoff [Å]      default: 0.7
  --elements Fe,O     Include only these elements default: all
  --last-n  INT       Use only the last N frames
  --ncore   INT       Parallel threads
  -o STEM             Output name stem            default: radius.cube
  --outdir DIR        Write the cubes here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Example:
  ferro map radius -i traj.dump --elements Li --radius 0.7
  ferro map radius -i traj.dump --elements Li --radius 1.0 --nx 100 --ny 100 --nz 100

Full documentation:  ferro doc map radius"#
    );
}

pub fn print_cube_sdf() {
    println!(
        r#"ferro map sdf — Cluster Spatial Distribution Function
  Identifies Qn-type clusters (connected components of network-former atoms
  linked by bridging ligands), aligns each cluster to a reference via
  Kabsch rotation, and accumulates per-atom-type 3D probability density maps.
  Clusters with identical atom-type composition are grouped into the same
  family. The first cluster encountered per family is used as the reference.
  Outputs one Gaussian cube file per atom type per family.

  Atom-type labels:
    Former (e.g. P):  P0 / P1 / P2 / P3  (individual Qn connectivity)
    Ligand  (e.g. O): Of (free), On (non-bridging), Ob (bridging)
    Modifier (e.g. Zn): element symbol

  Output files:  <stem>_<atom_type>.cube           (single family)
                 <stem>_fam<N>_<atom_type>.cube     (multiple families)

Parameters:
  --qn         INT    Target Qn cluster level (0/1/2/3)       default: 3
  --former     ELEM   Network-former element                   default: P
  --ligand     ELEM   Ligand (bridging) element                default: O
  --cutoff-fl  FLOAT  Former-ligand bond cutoff [Å]           default: 2.4
  --modifier   ELEM   Modifier element (optional, e.g. Zn)
  --cutoff-ml  FLOAT  Modifier-ligand cutoff [Å]              default: 2.8
  --grid-res   FLOAT  Voxel size [Å]                          default: 0.1
  --sigma      FLOAT  Gaussian broadening sigma [voxels]       default: 1.5
  --padding    FLOAT  Grid boundary padding [Å]               default: 3.0
  --rmsd-warn  FLOAT  RMSD warning threshold [Å]              default: 0.5
  --last-n     INT    Use only the last N frames
  --ncore      INT    Parallel threads
  -o STEM             Output stem (no extension)              default: sdf
  --outdir DIR        Write the cubes here (created if missing)
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Example:
  ferro map sdf -i traj.dump --qn 3
  ferro map sdf -i traj.dump --qn 2 --modifier Zn --cutoff-ml 2.8 -o q2_sdf
  ferro map sdf -i traj.dump --qn 1 --grid-res 0.05 --sigma 2.0 --last-n 500

Full documentation:  ferro doc map sdf"#
    );
}

pub fn print_cube_chg_sdf() {
    println!(
        r#"ferro map chg-sdf — Averaged Charge-Density Cluster SDF
  Reads multiple QE pp.x cube files (one per MD frame), identifies Qn
  clusters with the same logic as -m sdf, extracts a cubic sub-grid of
  the charge density centered on the cluster anchor, applies the Kabsch
  rotation to align the sub-grid to a common reference frame, and
  accumulates the averaged charge density.

  Input: --cubes <file1.cube> <file2.cube> ...
  Each cube file contains both atomic structure and charge density.
  All cube files must have the same grid resolution (i.e. same QE cutoff).

  Output values are in ChargeGrid convention (ρ_phys × V_cell).
  The output cube file can be visualised directly in VESTA or VMD.

Parameters:
  --cubes      FILE...  QE pp.x cube files (required, one per frame)
  --qn         INT      Target Qn cluster level (0/1/2/3)      default: 2
  --former     ELEM     Network-former element                  default: P
  --ligand     ELEM     Ligand (bridging) element               default: O
  --cutoff-fl  FLOAT    Former-ligand bond cutoff [Å]          default: 2.4
  --modifier   ELEM     Modifier element (optional, e.g. Zn)
  --cutoff-ml  FLOAT    Modifier-ligand cutoff [Å]             default: 2.8
  --chg-padding FLOAT   Sub-grid boundary margin [Å]           default: 6.0
  --rmsd-warn  FLOAT    RMSD warning threshold [Å]             default: 0.5
  --ncore      INT      Parallel threads
  -o STEM               Output stem (no extension)             default: chg_sdf
  --outdir DIR          Write the cubes here (created if missing)

Example:
  ferro map chg-sdf --cubes frame*.cube --qn 2 --former P --ligand O -o Q2_avg
  ferro map chg-sdf --cubes f1.cube f2.cube --qn 0 --chg-padding 5.0

Full documentation:  ferro doc map chg-sdf"#
    );
}

pub fn print_dataset_overview() {
    println!(
        r#"ferro dataset — Machine-learning training sets

  Three steps kept as separate commands, because the first one is expensive and
  its output is the copy you back up:

  collect    AIMD output -> DeePMD system directories        (implemented)
  filter     quality selection on an existing dataset        (implemented)
  merge      combine same-composition datasets, resize sets  (implemented)

Run a subcommand with no -i for its full page:
  ferro dataset collect"#
    );
}

pub fn print_dataset_collect() {
    println!(
        r#"ferro dataset collect — AIMD output -> DeePMD system directories

  Reads CP2K MD output (the stdout log, with coordinates, forces and stress all
  printed to __STD_OUT__) and writes one DeePMD system per input DIRECTORY.

Parameters:
  -i, --input  FILE...    AIMD output files; glob patterns allowed
  -o, --outdir DIR        Output root (required)
      --overwrite         Allow writing into an existing non-empty directory

Output layout:
  <outdir>/<dir below the shared ancestor>/
    type.raw  type_map.raw  set.000/coord|box|energy|force|virial .npy

  -i 'run*/*.out' -o sets            -> sets/run1/, sets/run2/
  -i '/s/a/md/x.out' '/s/b/md/x.out' -> sets/a/md/, sets/b/md/
  -i '*.out'       (one directory)   -> sets/ itself

One system per directory:
  The .out files of one directory are the restart segments of one run, so they
  are reassembled into ONE system. That is the line between the two commands:
  collect puts back together one run, merge combines different runs.

  Files are ordered by their first MD| Step number; overlapping frames are kept,
  and every source file's step span is printed so the overlap stays visible.
  Two compositions in one directory is an error, not a frame-dropping event.

Examples:
  ferro dataset collect -i 'run*/*.out' -o data
  ferro dataset collect -i md1.out md2.out -o /scratch/train

Full documentation:  ferro doc dataset collect"#
    );
}

pub fn print_dataset_filter() {
    println!(
        r#"ferro dataset filter — Drop low-quality frames from a dataset

  Reads DeePMD system directories, keeps the good frames, writes them as a new
  dataset. The input is never modified.

Parameters:
  -i, --input  DIR...     System directories, or a directory holding them
                          (searched recursively)
  -o, --outdir DIR        Output root; each system is rebuilt under its path
                          relative to -i. OMIT for a read-only run
  -f, --f-max  EV_PER_A   Largest force magnitude allowed             [20.0]
  -s, --s-max  GPA        Largest |stress component| allowed          [10.0]
      --start  N          First SURVIVING frame to take (0-based)        [0]
      --end    N          Last SURVIVING frame (0-based, INCLUSIVE)   [last]
      --stride N          Take every Nth surviving frame                 [1]
  -N, --number N          Take this many surviving frames, spread evenly
      --oo-min [DMIN]     Drop frames whose smallest O-O distance is below
                          this; bare --oo-min uses 2.0 A
      --al6 [RCUT]        Keep only frames holding a 6-coordinated Al; bare
                          --al6 takes the cutoff from the Al-O RDF
      --shuffle           Shuffle the kept frames, after every other step
      --seed   N          Seed for --shuffle                          [666]
      --set-size N        Frames per output set; 0 = one set          [400]
      --overwrite         Allow writing into an existing non-empty directory

The funnel:
  all -> |F|max -> |sigma|max -> min d(O-O) -> Al6 -> [start:end:stride|N] -> shuffle

  A threshold of 0 switches that criterion off — an explicit zero says "do not
  judge", which no small positive number can express.
  --start/--end/--stride/-N count SURVIVING frames, not original indices.

Output:
  <outdir>/<path relative to -i>/   the filtered systems
  <outdir>/filter_*.csv            funnel, criteria, overlap + 4 diagnostics
  Without -o nothing is written; every table is printed instead.

Examples:
  ferro dataset filter -i raw                       # look, write nothing
  ferro dataset filter -i raw -o clean -f 15 -s 8
  ferro dataset filter -i raw -o clean -N 500 --set-size 250 --shuffle

Full documentation:  ferro doc dataset filter"#
    );
}

pub fn print_dataset_merge() {
    println!(
        r#"ferro dataset merge — Combine datasets of the same composition

  Reads DeePMD system directories, groups them by what they actually contain,
  and writes one merged system per composition.

Parameters:
  -i, --input  DIR...     System directories to combine (globs allowed)
  -o, --outdir DIR        Output root; one directory per composition (required)
      --mode   MODE       shuffle | by-source                    [shuffle]
      --seed   N          Shuffle seed; by-source ignores it          [666]
      --set-size N        Frames per output set; 0 = one set          [400]
      --suffix EXT        Force this suffix; default inherits a shared one
      --overwrite         Allow writing into an existing non-empty directory

Grouping:
  NOT by directory name — `init.011` says nothing reliable about its contents.
  Systems group by their per-atom element sequence and are written as
  <natoms>_<formula>, e.g. 7_Al2O4Zn (subscripts are actual counts).

The two modes:
  shuffle     One composition concatenated, shuffled with --seed, cut into sets.
              Every set holds a mix of whatever conditions went in.
  by-source   No mixing, no shuffling; each system is cut into sets ON ITS OWN,
              so every set.NNN stays traceable to one condition. The mapping is
              written to sets_source.txt.

  In both modes the remainder is spread evenly rather than left at the end:
  500 frames at --set-size 400 give 250 + 250, not 400 + 100.

Examples:
  ferro dataset merge -i 'data/*.train' -o merged
  ferro dataset merge -i 'data/*.train' -o merged --mode by-source
  ferro dataset merge -i 'clean/*' -o merged --set-size 250 --suffix .train

Full documentation:  ferro doc dataset merge"#
    );
}
