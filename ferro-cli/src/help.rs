use crate::args::{corr::CorrMode, cube::CubeCliMode, traj::TrajMode};

pub fn print_traj_help(mode: &TrajMode) {
    match mode {
        TrajMode::Gr    => print_gr(),
        TrajMode::Sq    => print_sq(),
        TrajMode::Msd   => print_msd(),
        TrajMode::Angle => print_angle(),
    }
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
    }
}

// ─── mol-traj modes ──────────────────────────────────────────────────────────

fn print_gr() {
    println!(
        r#"mol-traj -m gr — Radial Distribution Function
  Computes g(r) for all atom-pair types, coordination number CN(r),
  and per-pair bond-length statistics (mean, std, count).
  Requires periodic cell (PBC) in the input file.

Parameters:
  --r-max  FLOAT   Max cutoff radius [Å]           default: 10.005
  --dr     FLOAT   Histogram bin width [Å]          default: 0.01
  --r-cut  FLOAT   CN integration cutoff [Å]        default: 2.3
  --last-n INT     Use only the last N frames
  --ncore  INT     Parallel threads (default: all cores)
  -o PATH          Output file                      default: gr.dat
                   Also writes: <stem>_cn.dat

Example:
  mol-traj -m gr -i traj.xyz
  mol-traj -m gr -i traj.dump --r-max 8.0 --r-cut 3.2 --last-n 500 -o result.dat"#
    );
}

fn print_sq() {
    println!(
        r#"mol-traj -m sq — Structure Factor S(q)
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

Example:
  mol-traj -m sq -i traj.xyz
  mol-traj -m sq -i traj.xyz --weighting xrd --q-max 20.0 -o sq_xrd.dat"#
    );
}

fn print_msd() {
    println!(
        r#"mol-traj -m msd — Mean Square Displacement
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
  mol-traj -m msd -i traj.xyz --dt 2.0
  mol-traj -m msd -i traj.dump --elements Li --dt 1.0 --last-n 2000"#
    );
}

fn print_angle() {
    println!(
        r#"mol-traj -m angle — Bond Angle Distribution
  Computes P(θ) for all A-B-C triplets within cutoff distances.
  B is the central atom; A and C are its neighbors.

Parameters:
  --r-cut-ab FLOAT  A-to-B distance cutoff [Å]   default: 2.3
  --r-cut-bc FLOAT  C-to-B distance cutoff [Å]   default: 2.3
  --d-angle  FLOAT  Histogram bin width [°]       default: 0.1
  --last-n   INT    Use only the last N frames
  --ncore    INT    Parallel threads
  -o PATH           Output file                   default: angle.dat

Example:
  mol-traj -m angle -i traj.xyz
  mol-traj -m angle -i traj.xyz --r-cut-ab 2.0 --r-cut-bc 2.0"#
    );
}

// ─── mol-corr modes ──────────────────────────────────────────────────────────

fn print_vacf() {
    println!(
        r#"mol-corr -m vacf — Velocity Autocorrelation Function
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
  mol-corr -m vacf -i traj.dump --dt 2.0
  mol-corr -m vacf -i traj.dump --elements O --last-n 1000"#
    );
}

fn print_rotcorr() {
    println!(
        r#"mol-corr -m rotcorr — Rotational Correlation Function
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
  mol-corr -m rotcorr -i traj.xyz --center O --neighbor H
  mol-corr -m rotcorr -i traj.dump --center O --neighbor H --dt 2.0"#
    );
}

fn print_vanhove() {
    println!(
        r#"mol-corr -m vanhove — Van Hove Self-Correlation Function
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
  mol-corr -m vanhove -i traj.xyz --tau 100
  mol-corr -m vanhove -i traj.dump --elements Li --tau 500 --dt 2.0"#
    );
}

// ─── mol-cube modes ──────────────────────────────────────────────────────────

fn print_cube_density() {
    println!(
        r#"mol-cube -m density — Spatial Number Density
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
  mol-cube -m density -i traj.dump
  mol-cube -m density -i traj.dump --nx 100 --ny 100 --nz 100 --elements Li"#
    );
}

fn print_cube_velocity() {
    println!(
        r#"mol-cube -m velocity — Spatial Velocity Distribution
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
  mol-cube -m velocity -i traj.dump --nx 80 --ny 80 --nz 80"#
    );
}

fn print_cube_force() {
    println!(
        r#"mol-cube -m force — Spatial Force Distribution
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
  mol-cube -m force -i traj.dump --elements O"#
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
