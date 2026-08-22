
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

/// `ferro convert` with no `-i`: the read/write format matrix.
pub fn print_convert() {
    println!(
        r#"ferro convert — Structure / trajectory format conversion
  Reads one file, writes another. Both formats come from the file names;
  there is no --from / --to flag.

Supported formats:
{}

Parameters:
  -i, --input  FILE       Input file  (format from its name)
  -o, --output FILE       Output file (format from its name) — a full PATH here,
                          unlike the analysis commands where -o is a suffix
      --start  N          First frame to take      (0-based, inclusive) [0]
      --end    N          Last frame to take       (0-based, INCLUSIVE) [last]
      --stride N          Take every Nth frame within [start, end]      [1]
      --number N          Take this many frames, spread evenly over
                          [start, end] with both ends included
      --metal-units       Read/write LAMMPS dump in metal units
                          (velocities Å/ps, forces eV/Å; default is real units)
  -h, --help              Parameter table (this page shows the formats)

Selecting frames:
  --end is INCLUSIVE and 0-based, matching the frame numbers `ferro info` prints:
  a run whose last frame `info` calls `Frame 4` is fully covered by `--end 4`.

  --stride and --number are two ways of saying the same thing and cannot be
  combined: --stride is a spacing, --number is a total. --number always takes
  both endpoints, because the last frame of a run is usually the most
  equilibrated one and a fixed-step walk would systematically drop it.
  Asking for more frames than exist gives every frame once, never duplicates.

  ferro convert -i traj.dump -o sub.extxyz --start 100          # drop equilibration
  ferro convert -i traj.dump -o sub.extxyz --start 100 --end 199
  ferro convert -i traj.dump -o POSCAR --stride 50              # every 50th frame
  ferro convert -i traj.dump -o conf.lmp --number 20            # 20 spread evenly

How many files come out:
  One, if the target format holds a trajectory (see the Frames column above).
  One PER FRAME, if it holds a single structure — writing 20 frames to POSCAR
  can only be 20 files. The frame number is inserted before the extension, and
  it is the index in the ORIGINAL trajectory, so products trace back to it:

    -o POSCAR   --stride 2   ->  POSCAR_0000     POSCAR_0002     POSCAR_0004
    -o conf.lmp --number 3   ->  conf_0000.lmp   conf_0002.lmp   conf_0004.lmp

  Zero-padded to at least 4 digits so `ls` sorts them in frame order. Selecting
  a single frame writes one file with no number, whatever the format.

What survives a conversion:
  Positions and cell always. Velocities, forces and energy only if BOTH sides
  carry them — converting a .dump to .xyz silently drops velocities, because
  plain XYZ has nowhere to put them. Use .extxyz to keep them.

  The element column is always written as a CLEAN element symbol, whatever
  `Atom::label` holds. Only `ferro net --export-traj` folds site labels back
  into the LAMMPS dump element column.

Example:
  ferro convert -i input.xyz -o output.pdb
  ferro convert -i input.cif -o POSCAR
  ferro convert -i CONTCAR -o final.cif
  ferro convert -i traj.lammpstrj -o traj.extxyz --metal-units"#,
        crate::io_dispatch::supported_formats()
    );
}

/// `ferro info` with no `-i`: what the summary reports.
pub fn print_info() {
    println!(
        r#"ferro info — Structure / trajectory summary
  Prints frame count, composition, cell parameters, volume and mass density.
  Reads every format `ferro convert` can read (run `ferro convert` for the list).

Parameters:
  -i, --input  FILE       Input file (format from its name)
      --metal-units       Read LAMMPS dump in metal units (velocities Å/ps,
                          forces eV/Å; default is real units)
  -h, --help              Parameter table

Reported per frame — the FIRST and the LAST only, not every frame:
  Atoms       total count + per-element composition
  Cell        a b c (Å) and α β γ (°); "none (non-periodic)" for molecules
  Volume      Å³
  Density     g/cm³ — Σ(atomic mass) / cell volume, from the element table
              unless the file carries explicit masses. Omitted entirely when
              there is no cell: density without a volume is undefined, and a
              placeholder would read like a measurement.
  PBC         per-axis periodicity flags
  Energy / Forces / Velocities   whether the frame carries them

  A trajectory whose volume drifts between the first and last frame is an NPT
  run; the density line will drift with it. For a mean ± σ over ALL frames,
  read the `# volume = <mean> +/- <std>` header any `ferro traj` product carries.

Unknown elements:
  Masses come from the built-in element table. A symbol that is not in it
  (a stray site label, or the "X" a short PDB line degrades to) falls back to
  1 amu, which drags the density DOWN without any other visible sign — so the
  density line is followed by a warning naming how many atoms fell back.
  Treat the number as invalid until that warning is gone.

Example:
  ferro info -i input.xyz
  ferro info -i traj.lammpstrj
  ferro info -i CHGCAR"#
    );
}

/// `ferro bader` with no `-i`: methods, outputs, and the file-name collision.
pub fn print_bader() {
    println!(
        r#"ferro bader — Bader charge partitioning from a DFT charge density
  Partitions the charge density into atomic basins along its zero-flux surfaces,
  and reports the charge, volume and surface distance of each.

Input (format from the file name):
  .cube                   Gaussian / QE pp.x cube
  anything else           VASP CHGCAR

Parameters:
  -i, --input  FILE       Charge density file
  -m, --method NAME       ongrid | neargrid | offgrid | weight   default: neargrid
  -r, --refine INT        Edge refinement: -1 = auto, -2 = single pass,
                          N = N passes                           default: -1
  -v, --vacval FLOAT      Vacuum density threshold [e/Å³]        default: 1e-3
  -h, --help              Parameter table

Choosing a method:
  neargrid   Default. Gradient ascent with the accumulated off-grid correction,
             then edge refinement. Accurate on ordinary cells.
  ongrid     Cheapest, steepest-ascent between grid points only. Basin surfaces
             are staircased, so charges are systematically off by a little.
  offgrid    Interpolated gradients — slower, no grid bias.
  weight     Yu-Trinkle: every grid point is SPLIT between basins by flux
             weight instead of assigned whole. Use it for strongly skewed
             (non-orthogonal) cells, where the gradient direction the on/near
             grid methods use carries a known approximation error.

Output — three Henkelman-format .dat files, named after the INPUT file stem:
  <stem>_ACF.dat          per atom: charge, volume, min distance to the surface
  <stem>_BCF.dat          per Bader volume: charge, volume, position
  <stem>_AVF.dat          atom -> Bader volume index

  These stay in the Henkelman layout (not csv) because external tools parse them.

  CAUTION: they are written to the CURRENT DIRECTORY and there is no --outdir
  yet. VASP names every charge density CHGCAR, so running two systems from the
  same working directory writes CHGCAR_ACF.dat twice — the second run silently
  overwrites the first. Until --outdir lands, cd into each system's directory
  (or rename the inputs) rather than running them side by side.

Example:
  ferro bader -i CHGCAR
  ferro bader -i charge.cube
  ferro bader -i CHGCAR --method weight
  ferro bader -i CHGCAR --method neargrid --refine 3 --vacval 1e-4"#
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

Batch input:
  -i takes several files and expands glob patterns itself — quote them:
    ferro traj gr -i 'runs/*/prod.dump' -a P -b O -o scan
  Each input is analysed on its own; results stack into ONE csv with a `file` column.
  A failed input is skipped, listed in the output's [inputs] block, and sets exit code 1.

Output naming:
  <outdir>/<command>[_<table>][_<label>]_<suffix>.csv
  --outdir DIR   where every product goes (created if missing; default: current dir)
  -o SUFFIX      batch tag, chosen by you
  <label>        what was analysed — filled in from the type selection, so
                 `traj gr -a P -b O` writes gr_P-O.csv and no selection writes
                 gr_all.csv. Label before suffix, so `ls gr_P-O_*` lists one pair
                 across every batch.

Help:
  Any command typed without -i prints its own page: what it computes, the
  parameters, and the shape of what it writes. `-h` gives the short parameter
  table instead."#,
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
  Computes g(r) and CN(r) for every ordered pair of types, into one file.
  Requires periodic cell (PBC) in the input file.

  g(r) is symmetric: A-B and B-A hold identical values.
  CN(r) is directed:  A-B is the average number of B around each A,
                      so A-B and B-A generally differ.

Selecting a pair (centre first, neighbour second):
  -a ELEM  -b ELEM        by element, e.g. -a P -b O  -> O around each P
  -x LABEL -y LABEL       by site label, e.g. -x P_2 -y O_b
  (the two groups are mutually exclusive; omit both to get every pair)

  Site labels come from the LAMMPS dump element column written as
  <Element>_<suffix> (P_2, O_b, Al_4 — what `ferro net --export-traj` writes);
  they are split into element + label on read, so -a/-b keep working on the
  plain element.

  NOTE: -x/-y hold over a SINGLE frame only. g(r) requires a fixed particle count
  per type, and labels shift as the run evolves, so a multi-frame labelled
  trajectory is rejected. Use --last-n 1, or select by element.

Parameters:
  --r-min  FLOAT          Min cutoff radius [Å]                  default: 0.001
  --r-max  FLOAT          Max cutoff radius [Å]                  default: 10.005
                          (clamped to half the smallest interplanar spacing)
  --dr     FLOAT          Histogram bin width [Å]                default: 0.002
  --last-n INT            Use only the last N frames
  --ncore  INT            Parallel threads (default: all cores)
  -o SUFFIX               Batch tag  -> gr_<pair>_<suffix>.csv
  --outdir DIR            Write products here (created if missing)
  --plot                  PNG next to the data file (needs a pair)

File name — gr_<pair>[_<suffix>].csv, the pair as you wrote it:
  -a P -b O  -> gr_P-O.csv        -a O -b P  -> gr_O-P.csv   (different data: cn
                                                              is directed)
  no pair    -> gr_all.csv        -x P_3 -y O_b -> gr_P_3-O_b.csv

Output — long format, one row per (file, r, pair):
  file  r  center  neighbor  gr  cn

  The pair lives in data columns, not column names, so runs with different element
  sets stack without any column alignment. Omitting -a/-b just adds rows (every
  ordered pair), never columns. gr is symmetric; cn is directed (center -> neighbor).

Example:
  ferro traj gr -i traj.dump
  ferro traj gr -i traj.dump -a P -b O --r-max 8.0 --plot
  ferro traj gr -i 'runs/*/prod.dump' -a P -b O -o scan     # 一次跑一批
  ferro traj gr -i traj.dump -a P -b O --outdir results/700K
  ferro traj gr -i traj.dump -x P_2 -y O_b --last-n 1"#
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

No type selection: EVERY pair is written, always.
  -a/-b and -x/-y were removed. The primary product is the pair of totals; the
  partials are a decomposition that sums back to them, and keeping one pair hides
  exactly that closure. Filter columns in pandas instead — the file is small.
  Label-resolved partials went with them: a site label rarely carries enough atoms
  for its partial to show a signal.

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
  -o SUFFIX           Batch tag -> sq_<suffix>.csv   default: sq.csv
  --outdir DIR        Write products here (created if missing)
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
  ferro traj sq -i traj.dump --outdir results/700K"#
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

File name — msd_<elements>[_<suffix>].csv, elements sorted:
  --elements P,O -> msd_O-P.csv     (so does --elements O,P: it is a set, and the
                                     same data must not land under two names)
  no --elements  -> msd_all.csv

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
  -o SUFFIX                     Batch tag -> angle_<triplet>_<suffix>.csv
  --outdir DIR                  Write products here (created if missing)
  --plot                        Generate PNG and open in viewer

File name — angle_<triplet>[_<suffix>].csv, the triplet as you wrote it:
  -a O -b P -c O -> angle_O-P-O.csv        no triplet -> angle_all.csv

Example:
  ferro traj angle -i traj.dump
  ferro traj angle -i traj.dump -a O -b P -c O --r-cut-ab 2.0 --r-cut-bc 2.0
  ferro traj angle -i traj.dump -x O_b -y P_2 -z O_n --last-n 1
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
  -o SUFFIX             Batch tag -> vacf_<elements>_<suffix>.csv
  --outdir DIR          Write products here (created if missing)

File name — vacf_<elements>[_<suffix>].csv, elements sorted; vacf_all.csv without
  --elements.

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
  -o SUFFIX           Batch tag -> rotcorr_<centre>-<neighbour>_<suffix>.csv
  --outdir DIR        Write products here (created if missing)

File name — rotcorr_<centre>-<neighbour>[_<suffix>].csv; both are required, so this
  one never falls back to "all". --center O --neighbor H -> rotcorr_O-H.csv

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
  -o SUFFIX             Batch tag -> vanhove_<elements>_<suffix>.csv
  --outdir DIR          Write products here (created if missing)

File name — vanhove_<elements>[_<suffix>].csv, elements sorted; vanhove_all.csv
  without --elements.

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
  -o STEM             Output name stem            default: density.cube
  --outdir DIR        Write the cubes here (created if missing)

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
  -o STEM             Output name stem            default: velocity.cube
  --outdir DIR        Write the cubes here (created if missing)

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
  -o STEM             Output name stem            default: force.cube
  --outdir DIR        Write the cubes here (created if missing)

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
  -o STEM             Output name stem            default: radius.cube
  --outdir DIR        Write the cubes here (created if missing)

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
  -o STEM             Output stem (no extension)              default: sdf
  --outdir DIR        Write the cubes here (created if missing)

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
  -o STEM               Output stem (no extension)             default: chg_sdf
  --outdir DIR          Write the cubes here (created if missing)

Example:
  ferro map chg-sdf --cubes frame*.cube --qn 2 --former P --ligand O -o Q2_avg
  ferro map chg-sdf --cubes f1.cube f2.cube --qn 0 --chg-padding 5.0"#
    );
}
