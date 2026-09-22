
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
  -o, --output PATH   Output file, or a directory written with a trailing / to
                      hold the default name (job.gjf / job.inp / pw.in)
      --mkdir         Create -o's directory without asking
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

Output:
  One pw.x input, the path -o names (pw.in by default). Only frame 0 is used.

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
  --auto-spin         Guess multiplicity from the structure: magmom sum, then
                      oxidation state + Hund, then the electron-parity floor

Output:
  One .gjf, the path -o names (job.gjf by default). Only frame 0 is used.

Examples:
  ferro job -s gaussian -i mol.xyz
  ferro job -s gaussian -i mol.xyz -m PBE0 -b def2-TZVP -o sp.gjf
  ferro job -s gaussian -i FeCl3.xyz --auto-spin      # high-spin multiplicity
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

What CP2K prints:
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

Output:
  One .inp, the path -o names (job.inp by default). Only frame 0 is used. The
  per-element basis and pseudopotential names come from a 2829-entry database;
  --cp2k-basis only picks the family.

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
  -o, --output FILE       Output file; may include directories that do not exist
                          yet (--mkdir creates them). A trailing / is refused:
                          the target format comes from the file name
      --start  N          First frame to take      (0-based, inclusive) [0]
      --end    N          Last frame to take       (0-based, INCLUSIVE) [last]
      --stride N          Take every Nth frame within [start, end]      [1]
      --number N          Take this many frames, spread evenly, both ends kept
      --mkdir             Create -o's directory without asking
      --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces
                          eV/Å); default is real units
  -h, --help              Short parameter table (this page adds the formats)

Output:
  One file when the target format holds a trajectory (Frames column above),
  one file PER FRAME when it holds a single structure (POSCAR, data, QE):
  POSCAR_0000, conf_0002.vasp — the index is the ORIGINAL frame number

Examples:
  ferro convert -i input.cif -o POSCAR
  ferro convert -i traj.dump -o sub.extxyz --start 100 --end 199
  ferro convert -i traj.dump -o out/conf.vasp --number 20 --mkdir

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

Output — stdout, the FIRST and the LAST frame only:
  atoms + composition, cell a b c / alpha beta gamma, volume, density (g/cm³),
  per-axis PBC flags, and whether energy / forces / velocities are present.
  No cell means no density line rather than a placeholder.

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

Parameters:
  -i, --input  FILE       Charge density file; .cube is Gaussian / QE pp.x,
                          anything else is read as a VASP CHGCAR
  -m, --method NAME       ongrid (cheapest, staircased basins) | neargrid
                          (accurate on ordinary cells) | offgrid (interpolated,
                          slower, no grid bias) | weight (Yu-Trinkle, for
                          strongly skewed cells)                   [neargrid]
  -r, --refine INT        Edge refinement: -1 auto, -2 single pass, N passes [-1]
  -v, --vacval FLOAT      Vacuum density threshold [e/Å³]             [1e-3]
  -o, --output DIR        Where the reports go    [the input's own directory]
  -s, --suffix SUFFIX     Tag -> <stem>_ACF_<suffix>.dat
      --mkdir             Create -o without asking (needed with no terminal)
  -h, --help              Short parameter table

Output:
  <stem>_ACF.dat  <stem>_BCF.dat  <stem>_AVF.dat, in Henkelman's layout because
  external tools parse them. They land NEXT TO THE INPUT rather than in the
  current directory: VASP calls every charge density CHGCAR, so two runs would
  otherwise overwrite one another. -o collects them elsewhere, -s tags them.

Examples:
  ferro bader -i CHGCAR
  ferro bader -i CHGCAR --method weight
  ferro bader -i CHGCAR --method neargrid --refine 3 --vacval 1e-4
  ferro bader -i run1/CHGCAR -o reports -s weight --mkdir

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

Conventions:
  -i takes several files and expands globs itself — quote them. Each input is
  analysed on its own and the results stack into ONE csv carrying a `file`
  column; a failed input is skipped and the exit code is 1.
  -o is a DIRECTORY for the commands that write several products, the output
  FILE for `convert` and `job`. -s tags a batch:
    <-o dir>/<command>[_<table>][_<label>]_<suffix>.csv
  A command typed without -i prints its own page; -h gives the short table.

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
  -o, --output DIR      Write every product here; --mkdir creates it unasked
  -s, --suffix SUFFIX   Batch tag -> <command>[_<table>][_<label>]_<suffix>.csv
      --last-n N        Use only the last N frames of the trajectory
      --ncore  N        Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units (velocities in A/ps)

Selecting types (gr / angle only):
  -a -b -c by element (-a P -b O), -x -y -z by site label (-x P_2 -y O_b); the
  two groups exclude each other. For gr the first is the centre, and the order
  matters for CN. sq takes no selection — every pair is always written."#
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
  -o, --output DIR      Write the cubes here; --mkdir creates it unasked
  -s, --suffix STEM     Output file stem (default depends on the command)
      --last-n N        Use only the last N frames
      --ncore  N        Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units (velocities in A/ps)

Output:
  One 3-D grid file per input, not a stacked table — nothing here stacks. The
  name carries the input stem (density_<stem>.cube) so inputs cannot overwrite
  one another. No summary table, no plot."#
    );
}

// ─── ferro traj: gr / sq / msd / angle ──────────────────────────────────────

pub fn print_gr() {
    println!(
        r#"ferro traj gr — Radial Distribution Function and Coordination Number

  g(r) and CN(r) for every ordered pair of types, into one file. Needs a
  periodic cell. g(r) is symmetric (A-B = B-A) but CN(r) is DIRECTED, so the
  two directions of a pair generally differ.

Parameters:
  -a ELEM  -b ELEM        Pair by element, CENTRE first: -a P -b O gives O
                          around each P; omit both pairs for every pair
  -x LABEL -y LABEL       Select by site label: -x P_2 -y O_b. SINGLE frame
                          only; excludes -a/-b. Use --last-n 1
  --r-min  FLOAT          Min cutoff radius [Å]                     [0.001]
  --r-max  FLOAT          Max cutoff radius [Å]                    [10.005]
                          (clamped to half the smallest interplanar spacing)
  --dr     FLOAT          Histogram bin width [Å]                   [0.002]
  --last-n INT            Use only the last N frames
  --ncore  INT            Parallel threads                    [all cores]
  -o DIR                  Output directory; --mkdir creates it unasked
  -s SUFFIX               Batch tag  -> gr_<pair>_<suffix>.csv
  --metal-units           LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)
  --plot                  PNG next to the data file (needs a pair)

Output:
  gr_<pair>[_<suffix>].csv, long format: file pair r g_r cn_r
  No selection -> gr_all.csv

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
  (Waasmaier-Kirfel) form factors or neutron scattering lengths. No type
  selection: EVERY pair is written, always. Filter columns in pandas.

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
  -o DIR                  Output directory; --mkdir creates it unasked
  -s SUFFIX               Batch tag -> sq_<suffix>.csv
  --metal-units           LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)
  --plot                  PNG next to the data file (weighted totals only)

Output:
  sq[_<suffix>].csv, WIDE format, one row per (file, q):
  file q total_xrd total_neutron, then three columns per pair
  (_sq unweighted, _xrd and _neutron carrying w_ij(q)*S_ij(q))

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

  MSD(t) = <|r(t₀+t) − r(t₀)|²> averaged over time origins, with the total and
  the per-axis (a/b/c) components.

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
  -o DIR                 Output directory; --mkdir creates it unasked
  -s SUFFIX              Batch tag -> msd_<elements>_<suffix>.csv
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  msd_<elements>[_<suffix>].csv, elements SORTED so one set has one name
  (--elements P,O and O,P both give msd_O-P.csv); no --elements -> msd_all.csv

Examples:
  ferro traj msd -i traj.xyz --dt 2.0
  ferro traj msd -i traj.dump --elements Li --dt 1.0 --last-n 2000
  ferro traj msd -i traj.dump --dt 1.0 --fit-range 0.3,0.8 --plot

Full documentation:  ferro doc traj msd"#
    );
}

pub fn print_angle() {
    println!(
        r#"ferro traj angle — Bond Angle Distribution

  P(θ) for all A-B-C triplets within cutoff distances. B is the central atom,
  and each geometric angle is counted ONCE: a PO4 tetrahedron gives 6 O-P-O
  angles, not 12.

Parameters:
  -a ELEM -b ELEM -c ELEM  Triplet by element,    B is the centre (all three
                           are required; excludes -x/-y/-z)
  -x LBL  -y LBL  -z LBL   Triplet by site label, Y is the centre
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
  -o DIR                  Output directory; --mkdir creates it unasked
  -s SUFFIX               Batch tag -> angle_<triplet>_<suffix>.csv
  --metal-units           LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)
  --plot                  PNG next to the data file

Output:
  angle_<triplet>[_<suffix>].csv, long format: file triplet theta count p
  The triplet is spelled as you wrote it (-a O -b P -c O -> angle_O-P-O.csv);
  no triplet -> angle_all.csv

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

  C_v(t) = <v(t₀)·v(t₀+t)> / <v²(t₀)> averaged over origins, plus its running
  integral (Green-Kubo D). Needs velocities in the input file.

Parameters:
  --dt       FLOAT      Timestep [fs]                 default: 1.0
  --shift    INT        Time-origin stride             default: 1
  --elements Fe,O,...   Include only these elements    default: all
  --last-n   INT        Use only the last N frames
  --tau      INT        Lag time in frames             default: half the run
  --ncore    INT        Parallel threads               default: all cores
  -o DIR                Output directory; --mkdir creates it unasked
  -s SUFFIX             Batch tag -> vacf_<elements>_<suffix>.csv
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  vacf_<elements>[_<suffix>].csv, elements sorted; vacf_all.csv without
  --elements

Examples:
  ferro traj vacf -i traj.dump --dt 2.0
  ferro traj vacf -i traj.dump --elements O --last-n 1000

Full documentation:  ferro doc traj vacf"#
    );
}

pub fn print_rotcorr() {
    println!(
        r#"ferro traj rotcorr — Rotational Correlation Function

  C₂(t) = <P₂(û(t₀)·û(t₀+t))> for molecular bond vectors. --center and
  --neighbor are required: they define the bond direction.

Parameters:
  --center    ELEM    Central atom element (required)   e.g. O
  --neighbor  ELEM    Neighbor atom element (required)  e.g. H
  --r-cut     FLOAT   Bond search cutoff [Å]            default: 1.2
  --dt        FLOAT   Timestep [fs]                     default: 1.0
  --shift     INT     Time-origin stride                default: 1
  --last-n    INT     Use only the last N frames
  --tau       INT     Lag time in frames                default: half the run
  --ncore     INT     Parallel threads                  default: all cores
  -o DIR              Output directory; --mkdir creates it unasked
  -s SUFFIX           Batch tag -> rotcorr_<centre>-<neighbour>_<suffix>.csv
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  rotcorr_<centre>-<neighbour>[_<suffix>].csv — never falls back to "all",
  since both ends are required (--center O --neighbor H -> rotcorr_O-H.csv)

Examples:
  ferro traj rotcorr -i traj.xyz --center O --neighbor H
  ferro traj rotcorr -i traj.dump --center O --neighbor H --dt 2.0

Full documentation:  ferro doc traj rotcorr"#
    );
}

pub fn print_vanhove() {
    println!(
        r#"ferro traj vanhove — Van Hove Self-Correlation Function

  Gs(r, τ), the distribution of atomic displacements over a fixed time lag τ.

Parameters:
  --tau      INT        Lag time in frames              default: half trajectory
  --dt       FLOAT      Timestep [fs]                  default: 1.0
  --shift    INT        Time-origin stride              default: 1
  --r-max    FLOAT      Max displacement [Å]           default: 10.0
  --dr       FLOAT      Bin width [Å]                  default: 0.01
  --elements Fe,O,...   Track only these elements       default: all
  --last-n   INT        Use only the last N frames
  --ncore    INT        Parallel threads                default: all cores
  -o DIR                Output directory; --mkdir creates it unasked
  -s SUFFIX             Batch tag -> vanhove_<elements>_<suffix>.csv
  --metal-units         LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  vanhove_<elements>[_<suffix>].csv, elements sorted; vanhove_all.csv without
  --elements

Examples:
  ferro traj vanhove -i traj.xyz --tau 100
  ferro traj vanhove -i traj.dump --elements Li --tau 500 --dt 2.0

Full documentation:  ferro doc traj vanhove"#
    );
}

// ─── ferro map ──────────────────────────────────────────────────────────────

pub fn print_cube_density() {
    println!(
        r#"ferro map density — Spatial Number Density

  Splits the box into nx*ny*nz voxels and averages the atom number density
  [atoms/Å³] over time.

Parameters:
  --nx INT            Grid points along a axis    default: 50
  --ny INT            Grid points along b axis    default: 50
  --nz INT            Grid points along c axis    default: 50
  --elements Fe,O     Count only these elements   default: all
  --last-n   INT      Use only the last N frames
  --ncore    INT      Parallel threads
  -o DIR              Output directory; --mkdir creates it unasked
  -s STEM             Output name stem            default: density.cube
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  <stem>.cube, one per input — a Gaussian cube VESTA and VMD read directly

Examples:
  ferro map density -i traj.dump
  ferro map density -i traj.dump --nx 100 --ny 100 --nz 100 --elements Li

Full documentation:  ferro doc map density"#
    );
}

pub fn print_cube_velocity() {
    println!(
        r#"ferro map velocity — Spatial Velocity Distribution

  Time-averaged speed |v| per voxel [Å/fs]. Needs velocities in the input.

Parameters:
  --nx INT            Grid points along a axis    default: 50
  --ny INT            Grid points along b axis    default: 50
  --nz INT            Grid points along c axis    default: 50
  --elements Fe,O     Include only these elements default: all
  --last-n   INT      Use only the last N frames
  --ncore    INT      Parallel threads
  -o DIR              Output directory; --mkdir creates it unasked
  -s STEM             Output name stem            default: velocity.cube
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  <stem>.cube, one per input

Examples:
  ferro map velocity -i traj.dump --nx 80 --ny 80 --nz 80

Full documentation:  ferro doc map velocity"#
    );
}

pub fn print_cube_force() {
    println!(
        r#"ferro map force — Spatial Force Distribution

  Time-averaged force magnitude |f| per voxel [eV/Å]. Needs forces in the input.

Parameters:
  --nx INT            Grid points along a axis    default: 50
  --ny INT            Grid points along b axis    default: 50
  --nz INT            Grid points along c axis    default: 50
  --elements Fe,O     Include only these elements default: all
  --last-n   INT      Use only the last N frames
  --ncore    INT      Parallel threads
  -o DIR              Output directory; --mkdir creates it unasked
  -s STEM             Output name stem            default: force.cube
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  <stem>.cube, one per input

Examples:
  ferro map force -i traj.dump --elements O

Full documentation:  ferro doc map force"#
    );
}

pub fn print_cube_radius() {
    println!(
        r#"ferro map radius — Hard-Sphere Spatial Occupancy Map

  Counts, per voxel, the (frame, atom) pairs with an atom within --radius of
  the voxel centre, under the minimum-image convention. The criterion is hard
  and binary, unlike the broadened count `map density` makes.

Parameters:
  --nx      INT       Grid points along a axis    default: 50
  --ny      INT       Grid points along b axis    default: 50
  --nz      INT       Grid points along c axis    default: 50
  --radius  FLOAT     Hard-sphere cutoff [Å]      default: 0.7
  --elements Fe,O     Include only these elements default: all
  --last-n  INT       Use only the last N frames
  --ncore   INT       Parallel threads
  -o DIR              Output directory; --mkdir creates it unasked
  -s STEM             Output name stem            default: radius.cube
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  <stem>.cube, one per input

Examples:
  ferro map radius -i traj.dump --elements Li --radius 0.7
  ferro map radius -i traj.dump --elements Li --radius 1.0 --nx 100 --ny 100 --nz 100

Full documentation:  ferro doc map radius"#
    );
}

pub fn print_cube_sdf() {
    println!(
        r#"ferro map sdf — Cluster Spatial Distribution Function

  Finds Qn clusters, aligns each to a reference by Kabsch rotation and
  accumulates a 3D probability density per atom type. Clusters of identical
  composition form one family, referenced to the first one seen.

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
  -o DIR              Output directory; --mkdir creates it unasked
  -s STEM             Output stem (no extension)              default: sdf
  --metal-units       LAMMPS dump in metal units (velocities Å/ps, forces eV/Å)

Output:
  <stem>_<atom_type>.cube, or <stem>_fam<N>_<atom_type>.cube with more than
  one family. Types are P0..P3 (former Qn), Of/On/Ob (free / non-bridging /
  bridging ligand) and the element symbol for a modifier.

Examples:
  ferro map sdf -i traj.dump --qn 3
  ferro map sdf -i traj.dump --qn 2 --modifier Zn --cutoff-ml 2.8 -o q2_sdf
  ferro map sdf -i traj.dump --qn 1 --grid-res 0.05 --sigma 2.0 --last-n 500

Full documentation:  ferro doc map sdf"#
    );
}

pub fn print_cube_chg_sdf() {
    println!(
        r#"ferro map chg-sdf — Averaged Charge-Density Cluster SDF

  Finds Qn clusters as `map sdf` does, cuts the charge density around each
  cluster anchor, Kabsch-aligns the sub-grids and averages them. All the input
  cubes must share one grid resolution, i.e. one QE cutoff.

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
  -o DIR                Output directory; --mkdir creates it unasked
  -s STEM               Output stem (no extension)             default: chg_sdf

Output:
  <stem>.cube — ONE cube from all the inputs together, unlike the rest of
  `map`. Values follow the ChargeGrid convention (rho_phys * V_cell).

Examples:
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

  Reads CP2K or VASP output and writes one DeePMD system per input DIRECTORY —
  the files of one directory are the restart segments of one run. The format
  comes from each file's own banner; --format overrides that.

Parameters:
  -i, --input  FILE...    AIMD output files; glob patterns allowed
  -o, --output DIR        Output root; default is beside each input
      --format FMT        Read the inputs as this format instead of trusting
                          the banner. cp2k/md | cp2k/sp | vasp/outcar |
                          vasp/xml                            [from banner]
      --type   WHAT       deepmd (DeePMD system) | inspect (diagnostics only)
                                                                    [deepmd]
      --mkdir             Create -o without asking (needed with no terminal)
      --overwrite         Allow writing into an existing non-empty directory

Output:
  <dir>.db/       one system per input directory: type.raw type_map.raw
                  set.NNN/coord|box|energy|force|virial .npy
  --type inspect  three diagnostic files in <dir>/ferro_inspect/ and NO
                  dataset, so -o is refused on that path

Examples:
  ferro dataset collect -i 'run*/*.out'
  ferro dataset collect -i 'run*/*.out' -o data
  ferro dataset collect -i 'md/*.out' --type inspect

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
  -o, --output DIR        Output root; each system is rebuilt under its path
                          relative to -i. OMIT for a read-only run
      --mkdir             Create -o without asking (needed with no terminal)
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
      --seed   N          Seed for --shuffle and for --ratio            [666]
      --set-size N        Frames per output set; 0 = one set          [400]
      --type   WHAT       deepmd | nep | extxyz                    [deepmd]
                          nep/extxyz write ONE .xyz per system instead of a
                          system directory; nep carries stress as virial=
      --ratio  A:B:C      Split train:valid:test, e.g. 8:1:1; two fields
                          mean train:test (9:1). Weights, not fractions
      --overwrite         Allow writing into an existing non-empty directory

Output:
  <-o>/<path below -i>.train/  filtered systems (one .xyz with --type nep)
  <-o>/filter_*.csv            funnel, criteria, overlap + 4 diagnostic tables
  Without -o nothing is written; every table is printed instead.

Examples:
  ferro dataset filter -i raw                       # look, write nothing
  ferro dataset filter -i raw -o clean -f 15 -s 8
  ferro dataset filter -i raw -o nep --type nep --ratio 8:1:1

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
  -o, --output DIR        Output root; one directory per composition (required)
      --mkdir             Create -o without asking (needed with no terminal)
      --mode   MODE       shuffle | by-source                    [shuffle]
      --seed   N          Shuffle seed; by-source ignores it          [666]
      --set-size N        Frames per output set; 0 = one set          [400]
      --type   WHAT       deepmd | nep | extxyz                    [deepmd]
                          nep/extxyz write ONE .xyz per group instead of a
                          system directory; nep carries stress as virial=
      --ratio  A:B:C      Split train:valid:test, e.g. 8:1:1; two fields
                          mean train:test (9:1). Weights, not fractions
      --suffix EXT        Force this suffix. Default: the one every input
                          shares, else .train; inputs from different parts of a
                          split are refused rather than silently relabelled
      --overwrite         Allow writing into an existing non-empty directory

Output:
  <-o>/<natoms>_<formula>/     one directory per composition, e.g. 7_Al2O4Zn
                               (or one .xyz with --type nep), suffixed .train
  --mode by-source also writes sets_source.txt; --ratio splits into
  .train / .valid / .test.

Examples:
  ferro dataset merge -i 'data/*.train' -o merged
  ferro dataset merge -i 'data/*.train' -o merged --mode by-source
  ferro dataset merge -i 'clean/*' -o nep --type nep --ratio 9:1

Full documentation:  ferro doc dataset merge"#
    );
}
