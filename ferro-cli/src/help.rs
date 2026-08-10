
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
      --metal-units   LAMMPS metal units for dump files"#
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

Parameters — task:
  --qe-task STR     Calculation type                    default: scf
    scf               Single-point SCF
    nscf              Non-self-consistent (after scf)
    bands             Band-structure run
    relax             Atomic relaxation (BFGS)
    vc-relax          Variable-cell relaxation
    md                Born-Oppenheimer MD
    vc-md             Variable-cell MD

Parameters — electronic structure:
  --qe-functional STR  DFT functional                   default: pbe
    pbe pbesol revpbe blyp scan r2scan pbe0 hse06
  --ecutwfc F       Plane-wave cutoff [Ry]               default: 50
  --smearing STR    Occupation smearing                  default: none
    none gaussian mp mv fd   (mp/mv recommended for metals)
  --kpoints K1 K2 K3  Monkhorst-Pack mesh (omit → Gamma)
  --pseudo-dir PATH Pseudopotential directory            default: ./pseudo

Charge / spin (shared):
  --charge INT        Override total charge
  --multiplicity INT  Override 2S+1 (→ nspin=2, tot_magnetization)
  --auto-spin         Guess spin from structure (default for qe);
                      nspin/tot_magnetization via guess_spin

Parameters — MD (--qe-task md|vc-md):
  --md-steps INT    Number of MD steps                   default: 10000
  --temperature F   Target temperature [K]               default: 298.15

Notes:
  ibrav = 0; cell from structure (CELL_PARAMETERS angstrom).
  Pseudopotentials referenced as <Element>.UPF in --pseudo-dir.

Examples:
  ferro job -s qe -i crystal.cif
  ferro job -s qe -i metal.cif --smearing mp --kpoints 8 8 8
  ferro job -s qe -i slab.xyz --qe-task relax --qe-functional scan
  ferro job -s qe -i Fe2O3.cif --auto-spin --kpoints 4 4 4 -o pw.in"#
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
  ferro job -s gaussian -i radical.xyz --charge 0 --multiplicity 2"#
    );
}

fn print_job_cp2k() {
    println!(
        r#"ferro job -s cp2k — CP2K input file (GPW/DFT, periodic systems)

Parameters — task:
  --task STR        Calculation type                    default: energy
    energy            Single-point energy
    force             Energy + forces
    geo-opt           Geometry optimisation (atoms)
    cell-opt          Geometry + cell optimisation
    md                Born-Oppenheimer molecular dynamics
    freq              Vibrational analysis

Parameters — electronic structure:
  --functional STR  DFT functional                      default: pbe
    pbe               GGA-PBE  (GTH-PBE pseudopotential)
    blyp              GGA-BLYP (GTH-BLYP pseudopotential)
    revpbe            GGA-revPBE
    pbesol            GGA-PBEsol
    pbe0              Hybrid PBE0  (25 % HF)
    b3lyp             Hybrid B3LYP (20 % HF, Gaussian definition)
    hse06             Range-separated HSE06
    scan              meta-GGA SCAN  (via LIBXC)
    r2scan            meta-GGA r²SCAN (via LIBXC)

  --cp2k-basis STR  Basis set                           default: dzvp-molopt-sr
    dzvp-molopt-sr    DZVP-MOLOPT-SR-GTH  (fast, good all-round)
    tzvp-molopt       TZVP-MOLOPT-GTH     (higher quality)
    tzv2p-molopt      TZV2P-MOLOPT-GTH    (highest quality MOLOPT)
    dzvp-gth          DZVP-GTH            (older GTH style)
    tzvp-gth          TZVP-GTH
    pob-dzvp          pob-DZVP            (all-electron, periodic)
    pob-tzvp          pob-TZVP            (all-electron, periodic)
  基组/赝势名经数据库按元素精确匹配（PBE/SCAN/全电子，含 q 价电子数）。

  --dispersion STR  Dispersion correction               default: none
    none  d3  d3bj

  --scf STR         SCF solver                          default: diag
    diag              Diagonalisation + Broyden  (metals, large systems)
    ot                Orbital Transform           (insulators / band-gap systems)

  --cutoff INT      Plane-wave cutoff [Ry]               default: 400
  --rel-cutoff INT  Relative cutoff [Ry]                 default: 50
  --smear           Enable Fermi-Dirac smearing (300 K)
  --pbc STR         Periodic boundary  xyz | z | none   (auto from cell)
  --kpoints K1 K2 K3  Monkhorst-Pack k-point mesh

Parameters — charge / spin (shared):
  --charge INT      Override total system charge
  --multiplicity INT  Override spin multiplicity 2S+1 (highest priority)
  --auto-spin       Guess multiplicity from structure:
                      magmom 求和 → 氧化态+Hund → 电子数奇偶下限
                      过渡金属按高自旋估计，结果需 DFT 验证

Parameters — output:
  --atom-charge STR Atomic charge scheme                 default: none
    none  mulliken  hirshfeld  hirshfeld-i
  --cube STR        Export cube file                     default: none
    none  density  elf  hartree
  --molden          Export Molden wavefunction file
  --project STR     CP2K project name                    default: ferro

Parameters — MD (--task md):
  --md-steps INT    Number of MD steps                   default: 10000
  --md-timestep F   Timestep [fs]                        default: 1.0
  --temperature F   Temperature [K]                      default: 298.15
  --thermostat STR  Thermostat                           default: csvr
    csvr              Canonical sampling (robust default)
    nose              Nosé-Hoover chain
    langevin          Langevin stochastic thermostat
    none              NVE (no thermostat)
  --traj-freq INT   Write trajectory every N steps       default: 100
  --barostat        Enable NPT barostat (flexible cell)

Examples:
  # Single-point PBE on a periodic glass structure
  ferro job -s cp2k -i glass.xyz

  # Geometry optimisation with DFT-D3(BJ)
  ferro job -s cp2k -i glass.xyz --task geo-opt --dispersion d3bj -o opt.inp

  # AIMD at 1500 K, NVT-CSVR, PBE-D3
  ferro job -s cp2k -i glass.xyz --task md --dispersion d3 \
         --temperature 1500 --md-steps 50000 --traj-freq 50 -o aimd.inp

  # Cell optimisation with PBE0 / OT / no dispersion
  ferro job -s cp2k -i crystal.cif --task cell-opt --functional pbe0 --scf ot

  # Mulliken charges + electron density cube
  ferro job -s cp2k -i mol.xyz --atom-charge mulliken --cube density

  # Auto-guess spin for a transition-metal oxide (high-spin estimate)
  ferro job -s cp2k -i Fe2O3.cif --auto-spin --smear"#
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
  net qn         Qn species / oxygen type / coordination number distributions
  net type       Per-atom type labels for one frame
  bader          Bader charge partitioning (ACF/BCF/AVF)

Structure I/O
  convert        Format conversion (xyz, pdb, cif, POSCAR, extxyz, lammps, qe, ...)
  info           Structure / trajectory information
  job            Quantum-chemistry input files (gaussian | cp2k | qe)

Batch input:
  -i takes several files and expands glob patterns itself — quote them:
    ferro traj gr -i 'runs/*/prod.dump' -a P -b O -o scan
  Each input is analysed on its own; results stack into ONE csv with a `file` column.
  A failed input is skipped, listed in the output's [inputs] block, and sets exit code 1."#,
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
  -o, --output SUFFIX   Output name suffix -> <command>[_<table>]_<suffix>.csv
      --last-n N        Use only the last N frames of the trajectory
      --ncore  N        Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units (velocities in A/ps)

Selecting types (gr / sq / angle):
  -a -b -c              by chemical element   (e.g. -a P -b O)
  -x -y -z              by site label         (e.g. -x P_0 -y O_b_P_P)
  The two groups are mutually exclusive. For gr and sq the first is the centre and the
  second the neighbour, and the order matters for CN."#
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

/// `ferro net` with no subcommand.
pub fn print_net_overview() {
    println!(
        r#"ferro net — Glass network topology

Usage:
  ferro net <COMMAND> -i <FILE> --<Former>-<Ligand>=<cutoff> [OPTIONS]

Commands:
  qn        Qn species / oxygen type / coordination number distributions (time averaged)
  type      Per-atom type labels for one frame; writes a structure file with -o

Cutoffs are given as flags naming the element pair, e.g. --P-O=2.3 --Zn-O=2.8.
Run `ferro net qn` without -i for the full option list and examples."#
    );
}

// ─── ferro traj: gr / sq / msd / angle ──────────────────────────────────────

pub fn print_gr() {
    println!(
        r#"ferro traj gr — Radial Distribution Function and Coordination Number
  Computes g(r) and CN(r) for every ordered pair of types, into one file.
  Requires periodic cell (PBC) in the input file.

  g(r) is symmetric: A-B and B-A hold identical values.
  CN(r) is directed:  A-B is the average number of B around each A,
                      so A-B and B-A generally differ.

Selecting a pair (centre first, neighbour second):
  -a ELEM  -b ELEM        by element, e.g. -a P -b O  -> O around each P
  -x LABEL -y LABEL       by site label, e.g. -x P_0 -y O_b_P_P
  (the two groups are mutually exclusive; omit both to get every pair)

  Site labels come from the LAMMPS dump element column written as
  <Element>_<suffix> (P_0, O_b_P_P, Zn_f); they are split into element +
  label on read, so -a/-b keep working on the plain element.

Parameters:
  --r-min  FLOAT          Min cutoff radius [Å]                  default: 0.001
  --r-max  FLOAT          Max cutoff radius [Å]                  default: 10.005
                          (clamped to half the smallest interplanar spacing)
  --dr     FLOAT          Histogram bin width [Å]                default: 0.002
  --last-n INT            Use only the last N frames
  --ncore  INT            Parallel threads (default: all cores)
  -o SUFFIX               Output suffix -> gr_<suffix>.csv        default: gr.csv
  --plot                  PNG next to the data file (needs a pair)

Output — long format, one row per (file, r, pair):
  file  r  center  neighbor  gr  cn

  The pair lives in data columns, not column names, so runs with different element
  sets stack without any column alignment. Omitting -a/-b just adds rows (every
  ordered pair), never columns. gr is symmetric; cn is directed (center -> neighbor).

Example:
  ferro traj gr -i traj.dump
  ferro traj gr -i traj.dump -a P -b O --r-max 8.0 --plot
  ferro traj gr -i 'runs/*/prod.dump' -a P -b O -o scan     # 一次跑一批
  ferro traj gr -i traj.dump -x P_0 -y O_b_P_P --last-n 500"#
    );
}

pub fn print_sq() {
    println!(
        r#"ferro traj sq — Structure Factor S(q)
  Computes S(q) via Fourier transform of g(r) (Faber-Ziman formalism),
  weighted by XRD (Waasmaier-Kirfel) form factors or neutron scattering lengths.

  Only the canonical half of the pairs is kept — S(q) is symmetric and has no
  directed counterpart, so B-A would duplicate A-B.

  Per pair three columns are written: _sq (unweighted), _xrd and _neutron
  (w_ij(q)*S_ij(q)). The weighted ones sum over pairs to total_xrd / total_neutron,
  so they decompose the experimentally comparable curve pair by pair.

Selecting a pair:
  -a ELEM  -b ELEM        by element
  -x LABEL -y LABEL       by site label — recommended here, since the pair count
                          grows quadratically with the number of labels and the
                          full table would run to hundreds of columns
  With a pair given only that pair's three columns are written, next to the totals.

Parameters:
  --q-min      FLOAT  Min q [Å⁻¹]                  default: 0.1
  --q-max      FLOAT  Max q [Å⁻¹]                  default: 25.0
  --dq         FLOAT  q bin width [Å⁻¹]            default: 0.02
  --weighting  ENUM   none | xrd | neutron | both   default: both
  --r-min      FLOAT  g(r) lower cutoff [Å]        default: 0.001
  --r-max      FLOAT  g(r) cutoff [Å]              default: 10.005
  --dr         FLOAT  g(r) bin width [Å]           default: 0.002
  --last-n     INT    Use only the last N frames
  --ncore      INT    Parallel threads (used in g(r) step)
  -o SUFFIX           Output suffix -> sq_<suffix>.csv   default: sq.csv
  --plot              PNG next to the data file (weighted totals only)

Output — wide format, one row per (file, q):
  file  q  total_xrd  total_neutron  then three columns per pair

  Wide, not long like gr: the primary product here is the pair of totals, which are
  one value per q. Inputs with different element sets contribute different pair
  columns; the gaps are left empty (NaN), never filled in.

Example:
  ferro traj sq -i traj.dump
  ferro traj sq -i traj.dump --weighting xrd --q-max 20.0 -o xrd
  ferro traj sq -i 'runs/*/prod.dump' -o scan
  ferro traj sq -i traj.dump -x P_0 -y O_b_P_P"#
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
  -o SUFFIX              Output suffix -> msd_<suffix>.csv   default: msd.csv

Example:
  ferro traj msd -i traj.xyz --dt 2.0
  ferro traj msd -i traj.dump --elements Li --dt 1.0 --last-n 2000
  ferro traj msd -i traj.dump --dt 1.0 --fit-range 0.3,0.8 --plot"#
    );
}

pub fn print_angle() {
    println!(
        r#"ferro traj angle — Bond Angle Distribution
  Computes P(θ) for all A-B-C triplets within cutoff distances.
  B is the central atom; A and C are its neighbours.

Selecting a triplet (all three required):
  -a ELEM  -b ELEM  -c ELEM      by element,    B is the centre
  -x LABEL -y LABEL -z LABEL     by site label, Y is the centre
  (the two groups are mutually exclusive)

Parameters:
  --r-cut-ab FLOAT              End-A-to-centre-B cutoff [Å]    default: 2.3
  --r-cut-bc FLOAT              End-C-to-centre-B cutoff [Å]    default: 2.3
                                A is the atom given by -a / -x, C the one by -c / -z.
                                Without a named triplet the two fall back to the
                                canonical (Z, symbol) order; when both ends are the
                                same type, min(ab, bc) applies to both.
  --angle-min FLOAT             Histogram lower edge [°]        default: 0.0
  --angle-max FLOAT             Histogram upper edge [°]        default: 180.0
                                Angles outside the window are discarded, not hidden.
  --d-angle  FLOAT              Histogram bin width [°]         default: 0.1
  --last-n   INT                Use only the last N frames
  --ncore    INT                Parallel threads
  -o SUFFIX                     Output suffix -> angle_<suffix>.csv  default: angle.csv
  --plot                        Generate PNG and open in viewer

Example:
  ferro traj angle -i traj.dump
  ferro traj angle -i traj.dump -a O -b P -c O --r-cut-ab 2.0 --r-cut-bc 2.0
  ferro traj angle -i traj.dump -x O_b_P_P -y P_0 -z O_f
  ferro traj angle -i traj.dump -a O -b P -c O --angle-min 90 --angle-max 130

Counting: each geometric angle once (a PO4 tetrahedron gives 6 O-P-O angles, not 12).
  code1/dump2analysis enumerates ordered end pairs, so its histogram is exactly twice
  this one when both ends are the same element; mean, std and peak positions agree.
  Its bins are also offset by half a bin — reproduce with --angle-min 0.05
  --angle-max 180.05 --d-angle 0.1."#
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
  -o SUFFIX             Output suffix -> vacf_<suffix>.csv   default: vacf.csv

Example:
  ferro traj vacf -i traj.dump --dt 2.0
  ferro traj vacf -i traj.dump --elements O --last-n 1000"#
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
  -o SUFFIX           Output suffix -> rotcorr_<suffix>.csv  default: rotcorr.csv

Example:
  ferro traj rotcorr -i traj.xyz --center O --neighbor H
  ferro traj rotcorr -i traj.dump --center O --neighbor H --dt 2.0"#
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
  -o SUFFIX             Output suffix -> vanhove_<suffix>.csv  default: vanhove.csv

Example:
  ferro traj vanhove -i traj.xyz --tau 100
  ferro traj vanhove -i traj.dump --elements Li --tau 500 --dt 2.0"#
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
  -o PATH             Output cube file            default: density.cube

Example:
  ferro map density -i traj.dump
  ferro map density -i traj.dump --nx 100 --ny 100 --nz 100 --elements Li"#
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
  -o PATH             Output cube file            default: velocity.cube

Example:
  ferro map velocity -i traj.dump --nx 80 --ny 80 --nz 80"#
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
  -o PATH             Output cube file            default: force.cube

Example:
  ferro map force -i traj.dump --elements O"#
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
  -o PATH             Output cube file            default: radius.cube

Example:
  ferro map radius -i traj.dump --elements Li --radius 0.7
  ferro map radius -i traj.dump --elements Li --radius 1.0 --nx 100 --ny 100 --nz 100"#
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
  -o PATH             Output stem (no extension)              default: sdf

Example:
  ferro map sdf -i traj.dump --qn 3
  ferro map sdf -i traj.dump --qn 2 --modifier Zn --cutoff-ml 2.8 -o q2_sdf
  ferro map sdf -i traj.dump --qn 1 --grid-res 0.05 --sigma 2.0 --last-n 500"#
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
  -o PATH               Output stem (no extension)             default: chg_sdf

Example:
  ferro map chg-sdf --cubes frame*.cube --qn 2 --former P --ligand O -o Q2_avg
  ferro map chg-sdf --cubes f1.cube f2.cube --qn 0 --chg-padding 5.0"#
    );
}
