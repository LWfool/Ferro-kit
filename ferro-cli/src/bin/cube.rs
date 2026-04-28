use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::cube::CubeCliMode,
    help::print_cube_help,
    io_dispatch::read_trajectory,
};
use ferro_analysis::{
    calc_cube_density, CubeDensityParams,
    calc_cube_radius, CubeRadiusParams,
    calc_cluster_sdf, ClusterSdfParams,
};
use ferro_io::{write_cube, LammpsUnits};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "fe-cube",
    about = "Spatial distribution maps: density | velocity | force | sdf  (run without -i for mode help)"
)]
struct Cli {
    /// Analysis mode
    #[arg(short = 'm', long, value_enum)]
    mode: CubeCliMode,

    /// Input trajectory file (omit to show mode-specific help)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output file stem (density/velocity/force) or stem prefix (sdf)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use only the last N frames
    #[arg(long)]
    last_n: Option<usize>,

    /// Parallel threads (default: all cores)
    #[arg(long)]
    ncore: Option<usize>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,

    // ── density / velocity / force / radius mode ─────────────────────────────

    /// Grid points along a axis  [density/velocity/force/radius]
    #[arg(long, default_value = "50")]
    nx: usize,

    /// Grid points along b axis  [density/velocity/force/radius]
    #[arg(long, default_value = "50")]
    ny: usize,

    /// Grid points along c axis  [density/velocity/force/radius]
    #[arg(long, default_value = "50")]
    nz: usize,

    /// Only include these elements, e.g. Fe,O  [density/velocity/force/radius]
    #[arg(long, value_delimiter = ',')]
    elements: Option<Vec<String>>,

    /// Hard-sphere radius cutoff [Å]  [radius]
    #[arg(long, default_value = "0.7")]
    radius: f64,

    // ── sdf mode ─────────────────────────────────────────────────────────────

    /// Target Qn cluster level (0/1/2/3)  [sdf]
    #[arg(long, default_value = "3")]
    qn: u8,

    /// Network-former element  [sdf]
    #[arg(long, default_value = "P")]
    former: String,

    /// Ligand (bridging) element  [sdf]
    #[arg(long, default_value = "O")]
    ligand: String,

    /// Former-ligand bond cutoff [Å]  [sdf]
    #[arg(long, default_value = "2.4")]
    cutoff_fl: f64,

    /// Modifier element (e.g. Zn); omit if none  [sdf]
    #[arg(long)]
    modifier: Option<String>,

    /// Modifier-ligand cutoff [Å]  [sdf]
    #[arg(long, default_value = "2.8")]
    cutoff_ml: f64,

    /// Voxel size [Å]  [sdf]
    #[arg(long, default_value = "0.1")]
    grid_res: f64,

    /// Gaussian broadening sigma [voxels]  [sdf]
    #[arg(long, default_value = "1.5")]
    sigma: f64,

    /// Grid boundary padding [Å]  [sdf]
    #[arg(long, default_value = "3.0")]
    padding: f64,

    /// RMSD warning threshold [Å]  [sdf]
    #[arg(long, default_value = "0.5")]
    rmsd_warn: f64,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    let input = match &args.input {
        Some(p) => p.clone(),
        None => {
            print_cube_help(&args.mode);
            return Ok(());
        }
    };

    if let Some(n) = args.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = args.last_n {
        traj = traj.tail(n);
    }

    match args.mode {
        CubeCliMode::Sdf    => run_sdf(&args, &traj),
        CubeCliMode::Radius => run_radius(&args, &traj),
        ref m               => run_density(m, &args, &traj),
    }
}

fn run_density(mode: &CubeCliMode, args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let default_name = format!("{}.cube", mode_name(mode));
    let out = args.output.as_deref().unwrap_or(Path::new(&default_name));

    let params = CubeDensityParams {
        nx: args.nx,
        ny: args.ny,
        nz: args.nz,
        elements: args.elements.clone(),
        mode: mode.clone().into(),
    };

    let result = calc_cube_density(traj, &params)
        .ok_or_else(|| anyhow!("Cube calc failed (missing cell, velocities, or forces?)"))?;

    write_cube(out.to_str().unwrap_or(&default_name), &result.cube)?;
    println!(
        "Cube ({}) -> {}  [{} frames, {} atoms]",
        mode_name(mode),
        out.display(),
        result.n_frames,
        result.n_atoms,
    );
    Ok(())
}

fn run_radius(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let out = args.output.as_deref().unwrap_or(std::path::Path::new("radius.cube"));

    let params = CubeRadiusParams {
        nx: args.nx,
        ny: args.ny,
        nz: args.nz,
        radius: args.radius,
        elements: args.elements.clone(),
    };

    let result = calc_cube_radius(traj, &params)
        .ok_or_else(|| anyhow!("Cube radius calc failed (missing cell?)"))?;

    write_cube(out.to_str().unwrap_or("radius.cube"), &result.cube)?;
    println!(
        "Cube (radius={:.3}Å) -> {}  [{} frames, {} atoms]",
        params.radius,
        out.display(),
        result.n_frames,
        result.n_atoms,
    );
    Ok(())
}

fn run_sdf(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let params = ClusterSdfParams {
        former: args.former.clone(),
        ligand: args.ligand.clone(),
        target_qn: args.qn,
        former_ligand_cutoff: args.cutoff_fl,
        modifier: args.modifier.clone(),
        modifier_cutoff: args.cutoff_ml,
        grid_res: args.grid_res,
        sigma: args.sigma,
        padding: args.padding,
        rmsd_warn_threshold: args.rmsd_warn,
    };

    let result = calc_cluster_sdf(traj, &params)
        .ok_or_else(|| anyhow!("No Q{} clusters found in trajectory", args.qn))?;

    let stem = args.output.as_deref()
        .and_then(|p| p.to_str())
        .unwrap_or("sdf");

    let multi_family = result.families.len() > 1;

    // Sort families by signature for deterministic file naming
    let mut families: Vec<_> = result.families.values().collect();
    families.sort_by(|a, b| a.signature.cmp(&b.signature));

    let mut total_files = 0usize;
    for (fam_idx, family) in families.iter().enumerate() {
        let fam_prefix = if multi_family {
            format!("{}_fam{}", stem, fam_idx)
        } else {
            stem.to_string()
        };

        // Sort atom types for deterministic output order
        let mut labels: Vec<_> = family.grids.keys().collect();
        labels.sort();

        for label in &labels {
            let cube = &family.grids[*label];
            let path = format!("{}_{}.cube", fam_prefix, label);
            write_cube(&path, cube)?;
            total_files += 1;
        }

        println!(
            "Family {:?}  ({} clusters, RMSD mean={:.3} max={:.3} Å, {} warnings)  → {}_*.cube",
            family.signature,
            family.n_clusters,
            family.rmsd_stats.mean,
            family.rmsd_stats.max,
            family.rmsd_stats.n_warned,
            fam_prefix,
        );
    }

    println!(
        "SDF Q{} done: {} frames, {} clusters total, {} cube files written",
        args.qn,
        result.n_frames,
        result.n_clusters_total,
        total_files,
    );
    Ok(())
}

fn mode_name(m: &CubeCliMode) -> &'static str {
    match m {
        CubeCliMode::Density  => "density",
        CubeCliMode::Velocity => "velocity",
        CubeCliMode::Force    => "force",
        CubeCliMode::Radius   => "radius",
        CubeCliMode::Sdf      => "sdf",
    }
}
