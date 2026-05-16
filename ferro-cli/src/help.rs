use crate::args::{corr::CorrMode, cube::CubeCliMode, traj::TrajMode};

pub fn print_fe_job_overview() {
    println!(
        r#"fe-job — Generate QC software input files

Usage:
  fe-job -s <SOFTWARE> -i <FILE> [OPTIONS]
  fe-job -s <SOFTWARE>              show software-specific parameters

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
        r#"fe-job -s qe — Quantum ESPRESSO pw.x input file

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
  fe-job -s qe -i crystal.cif
  fe-job -s qe -i metal.cif --smearing mp --kpoints 8 8 8
  fe-job -s qe -i slab.xyz --qe-task relax --qe-functional scan
  fe-job -s qe -i Fe2O3.cif --auto-spin --kpoints 4 4 4 -o pw.in"#
    );
}

fn print_job_gaussian() {
    println!(
        r#"fe-job -s gaussian — Gaussian 16/09 input file

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
  fe-job -s gaussian -i mol.xyz
  fe-job -s gaussian -i mol.xyz -m PBE0 -b def2-TZVP -o sp.gjf
  fe-job -s gaussian -i FeCl3.xyz --auto-spin            # 推断高自旋多重度
  fe-job -s gaussian -i radical.xyz --charge 0 --multiplicity 2"#
    );
}

fn print_job_cp2k() {
    println!(
        r#"fe-job -s cp2k — CP2K input file (GPW/DFT, periodic systems)

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
  fe-job -s cp2k -i glass.xyz

  # Geometry optimisation with DFT-D3(BJ)
  fe-job -s cp2k -i glass.xyz --task geo-opt --dispersion d3bj -o opt.inp

  # AIMD at 1500 K, NVT-CSVR, PBE-D3
  fe-job -s cp2k -i glass.xyz --task md --dispersion d3 \
         --temperature 1500 --md-steps 50000 --traj-freq 50 -o aimd.inp

  # Cell optimisation with PBE0 / OT / no dispersion
  fe-job -s cp2k -i crystal.cif --task cell-opt --functional pbe0 --scf ot

  # Mulliken charges + electron density cube
  fe-job -s cp2k -i mol.xyz --atom-charge mulliken --cube density

  # Auto-guess spin for a transition-metal oxide (high-spin estimate)
  fe-job -s cp2k -i Fe2O3.cif --auto-spin --smear"#
    );
}

pub fn print_fe_traj_overview() {
    println!(
        r#"fe-traj — Trajectory structural analysis

Usage:
  fe-traj -m <MODE> -i <FILE> [OPTIONS]
  fe-traj -m <MODE>              show mode-specific parameters

Modes:
  gr      Radial distribution function g(r) and coordination number CN(r)
  sq      Structure factor S(q) via Fourier transform of g(r)
  msd     Mean square displacement MSD(t), time-shift averaged
  angle   Bond angle distribution P(θ) for A-B-C triplets

Common options:
  -i, --input  PATH     Input trajectory file (xyz, dump, extxyz, pdb, …)
  -o, --output PATH     Output file (default name depends on mode)
      --last-n N        Use only the last N frames of the trajectory
      --ncore  N        Parallel threads (default: all cores)
      --plot            Generate a PNG plot and open it after calculation
      --metal-units     LAMMPS metal units (velocities in Å/ps)"#
    );
}

pub fn print_traj_help(mode: &TrajMode) {
    match mode {
        TrajMode::Gr    => print_gr(),
        TrajMode::Sq    => print_sq(),
        TrajMode::Msd   => print_msd(),
        TrajMode::Angle => print_angle(),
    }
}

pub fn print_fe_corr_overview() {
    println!(
        r#"fe-corr — Correlation functions

Usage:
  fe-corr -m <MODE> -i <FILE> [OPTIONS]
  fe-corr -m <MODE>              show mode-specific parameters

Modes:
  vacf      Velocity autocorrelation function C_v(t) + Green-Kubo diffusion
  rotcorr   Rotational correlation C₂(t) for molecular bond vectors
  vanhove   Van Hove self-correlation Gs(r, τ)

Common options:
  -i, --input  PATH     Input trajectory file (dump, xyz, extxyz, …)
  -o, --output PATH     Output file (default name depends on mode)
      --last-n N        Use only the last N frames
      --dt     FLOAT    Timestep between frames [fs]  (default: 1.0)
      --shift  INT      Time-origin stride            (default: 1)
      --metal-units     LAMMPS metal units (velocities in Å/ps)"#
    );
}

pub fn print_fe_cube_overview() {
    println!(
        r#"fe-cube — Spatial distribution maps

Usage:
  fe-cube -m <MODE> -i <FILE> [OPTIONS]
  fe-cube -m <MODE>              show mode-specific parameters

Modes:
  density   Time-averaged number density [atoms/Å³] per voxel
  velocity  Time-averaged speed |v| per voxel  (needs frame velocities)
  force     Time-averaged force magnitude |f| per voxel  (needs frame forces)
  radius    Hard-sphere spatial occupancy map
  sdf       Cluster spatial distribution function (Qn-type, Kabsch alignment)

Common options:
  -i, --input  PATH     Input trajectory file (dump, xyz, extxyz, …)
  -o, --output PATH     Output cube file / stem (default depends on mode)
      --last-n N        Use only the last N frames
      --ncore  N        Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units (velocities in Å/ps)"#
    );
}

pub fn print_corr_help(mode: &CorrMode) {
    match mode {
        CorrMode::Vacf    => print_vacf(),
        CorrMode::Rotcorr => print_rotcorr(),
        CorrMode::Vanhove => print_vanhove(),
    }
}

pub fn print_cube_help(mode: &CubeCliMode) {
    match mode {
        CubeCliMode::Density  => print_cube_density(),
        CubeCliMode::Velocity => print_cube_velocity(),
        CubeCliMode::Force    => print_cube_force(),
        CubeCliMode::Radius   => print_cube_radius(),
        CubeCliMode::Sdf      => print_cube_sdf(),
        CubeCliMode::ChgSdf   => print_cube_chg_sdf(),
    }
}

// ─── fe-traj modes ──────────────────────────────────────────────────────────

fn print_gr() {
    println!(
        r#"fe-traj -m gr — Radial Distribution Function
  Computes g(r) for all atom-pair types, coordination number CN(r),
  and per-pair bond-length statistics (mean ± std, count).
  Requires periodic cell (PBC) in the input file.

Parameters:
  --r-max  FLOAT          Max cutoff radius [Å]                  default: 10.005
  --dr     FLOAT          Histogram bin width [Å]                default: 0.01
  --r-cut  FLOAT          First-shell cutoff for pair stats [Å]  default: 2.3
  -a ELEM, -b ELEM        Show only one pair A-B (both required)
  --last-n INT            Use only the last N frames
  --ncore  INT            Parallel threads (default: all cores)
  -o PATH                 Output file                            default: gr.dat
                          Also writes: <stem>_cn.dat
  --plot                  Generate PNG and open in viewer

Example:
  fe-traj -m gr -i traj.xyz
  fe-traj -m gr -i traj.dump -a O -b P --r-max 8.0 --r-cut 2.0 --last-n 500"#
    );
}

fn print_sq() {
    println!(
        r#"fe-traj -m sq — Structure Factor S(q)
  Computes S(q) via Fourier transform of g(r) (Faber-Ziman formalism).
  Optionally applies XRD (Waasmaier-Kirfel) or neutron scattering weights.

Parameters:
  --q-max      FLOAT  Max q [Å⁻¹]                  default: 25.0
  --dq         FLOAT  q bin width [Å⁻¹]            default: 0.05
  --weighting  ENUM   none | xrd | neutron | both   default: both
  --r-max      FLOAT  g(r) cutoff [Å]              default: 10.005
  --dr         FLOAT  g(r) bin width [Å]           default: 0.01
  --last-n     INT    Use only the last N frames
  --ncore      INT    Parallel threads (used in g(r) step)
  -o PATH             Output file                   default: sq.dat
  --plot              Generate PNG and open in viewer

Example:
  fe-traj -m sq -i traj.xyz
  fe-traj -m sq -i traj.xyz --weighting xrd --q-max 20.0 -o sq_xrd.dat"#
    );
}

fn print_msd() {
    println!(
        r#"fe-traj -m msd — Mean Square Displacement
  Computes MSD(t) = <|r(t₀+t) − r(t₀)|²> averaged over time origins.
  Outputs total MSD and per-axis (a/b/c) components.

Parameters:
  --dt       FLOAT      Timestep between frames [fs]   default: 1.0
  --shift    INT        Time-origin stride             default: 1
  --elements Fe,O,...   Track only these elements      default: all
  --last-n   INT        Use only the last N frames
  --ncore    INT        Parallel threads
  -o PATH               Output file                    default: msd.dat

Example:
  fe-traj -m msd -i traj.xyz --dt 2.0
  fe-traj -m msd -i traj.dump --elements Li --dt 1.0 --last-n 2000"#
    );
}

fn print_angle() {
    println!(
        r#"fe-traj -m angle — Bond Angle Distribution
  Computes P(θ) for all A-B-C triplets within cutoff distances.
  B is the central atom; A and C are its neighbors.

Parameters:
  --r-cut-ab FLOAT              A-to-center-B bond cutoff [Å]   default: 2.3
  --r-cut-bc FLOAT              Center-B-to-C bond cutoff [Å]   default: 2.3
  --d-angle  FLOAT              Histogram bin width [°]         default: 0.1
  -a ELEM, -b ELEM, -c ELEM     Show only triplet A-B-C (all three required;
                                B is the center atom)
  --last-n   INT                Use only the last N frames
  --ncore    INT                Parallel threads
  -o PATH                       Output file                     default: angle.dat
  --plot                        Generate PNG and open in viewer

Example:
  fe-traj -m angle -i traj.xyz
  fe-traj -m angle -i traj.xyz -a O -b P -c O --r-cut-ab 2.0 --r-cut-bc 2.0"#
    );
}

// ─── fe-corr modes ──────────────────────────────────────────────────────────

fn print_vacf() {
    println!(
        r#"fe-corr -m vacf — Velocity Autocorrelation Function
  Computes C_v(t) = <v(t₀)·v(t₀+t)> / <v²(t₀)>, averaged over origins.
  Also outputs running integral (Green-Kubo diffusion coefficient).
  Requires frame.velocities in the input file.

Parameters:
  --dt       FLOAT      Timestep [fs]                 default: 1.0
  --shift    INT        Time-origin stride             default: 1
  --elements Fe,O,...   Include only these elements    default: all
  --last-n   INT        Use only the last N frames
  -o PATH               Output file                   default: vacf.dat

Example:
  fe-corr -m vacf -i traj.dump --dt 2.0
  fe-corr -m vacf -i traj.dump --elements O --last-n 1000"#
    );
}

fn print_rotcorr() {
    println!(
        r#"fe-corr -m rotcorr — Rotational Correlation Function
  Computes C₂(t) = <P₂(û(t₀)·û(t₀+t))> for molecular bond vectors.
  --center and --neighbor are required to define the bond direction.

Parameters:
  --center    ELEM    Central atom element (required)   e.g. O
  --neighbor  ELEM    Neighbor atom element (required)  e.g. H
  --r-cut     FLOAT   Bond search cutoff [Å]            default: 1.2
  --dt        FLOAT   Timestep [fs]                     default: 1.0
  --shift     INT     Time-origin stride                default: 1
  --last-n    INT     Use only the last N frames
  -o PATH             Output file                       default: rotcorr.dat

Example:
  fe-corr -m rotcorr -i traj.xyz --center O --neighbor H
  fe-corr -m rotcorr -i traj.dump --center O --neighbor H --dt 2.0"#
    );
}

fn print_vanhove() {
    println!(
        r#"fe-corr -m vanhove — Van Hove Self-Correlation Function
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
  -o PATH               Output file                    default: vanhove.dat

Example:
  fe-corr -m vanhove -i traj.xyz --tau 100
  fe-corr -m vanhove -i traj.dump --elements Li --tau 500 --dt 2.0"#
    );
}

// ─── fe-cube modes ──────────────────────────────────────────────────────────

fn print_cube_density() {
    println!(
        r#"fe-cube -m density — Spatial Number Density
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
  fe-cube -m density -i traj.dump
  fe-cube -m density -i traj.dump --nx 100 --ny 100 --nz 100 --elements Li"#
    );
}

fn print_cube_velocity() {
    println!(
        r#"fe-cube -m velocity — Spatial Velocity Distribution
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
  fe-cube -m velocity -i traj.dump --nx 80 --ny 80 --nz 80"#
    );
}

fn print_cube_force() {
    println!(
        r#"fe-cube -m force — Spatial Force Distribution
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
  fe-cube -m force -i traj.dump --elements O"#
    );
}

fn print_cube_radius() {
    println!(
        r#"fe-cube -m radius — Hard-Sphere Spatial Occupancy Map
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
  fe-cube -m radius -i traj.dump --elements Li --radius 0.7
  fe-cube -m radius -i traj.dump --elements Li --radius 1.0 --nx 100 --ny 100 --nz 100"#
    );
}

fn print_cube_sdf() {
    println!(
        r#"fe-cube -m sdf — Cluster Spatial Distribution Function
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
  fe-cube -m sdf -i traj.dump --qn 3
  fe-cube -m sdf -i traj.dump --qn 2 --modifier Zn --cutoff-ml 2.8 -o q2_sdf
  fe-cube -m sdf -i traj.dump --qn 1 --grid-res 0.05 --sigma 2.0 --last-n 500"#
    );
}

fn print_cube_chg_sdf() {
    println!(
        r#"fe-cube -m chg_sdf — Averaged Charge-Density Cluster SDF
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
  fe-cube -m chg_sdf --cubes frame*.cube --qn 2 --former P --ligand O -o Q2_avg
  fe-cube -m chg_sdf --cubes f1.cube f2.cube --qn 0 --chg-padding 5.0"#
    );
}
