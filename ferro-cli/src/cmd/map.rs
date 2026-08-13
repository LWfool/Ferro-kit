//! `ferro map` — spatial distribution maps.
//!
//! Grouped separately from `ferro traj` because the product is different: **one 3-D
//! `.cube` grid file per input**, not a table. That makes this the single exception to
//! the batch output rules — no stacked CSV, no `[inputs]` block, no plot — and the file
//! name must carry the input stem or a second input would overwrite the first.
//!
//! The traversal and failure handling are still the shared ones; error handling should
//! not exist in two versions.

use anyhow::{anyhow, Result};
use clap::{Args, Subcommand};
use ferro_analysis::{
    calc_chg_sdf, calc_cluster_sdf, calc_cube_density, calc_cube_radius, ChgSdfParams,
    ClusterSdfParams, CubeDensityParams, CubeRadiusParams,
};
use ferro_core::Trajectory;
use ferro_io::{read_cube_as_chg, write_cube};
use std::path::{Path, PathBuf};

use crate::args::common::CommonArgs;
use crate::args::cube::CubeCliMode;
use crate::batch;
use crate::help;

#[derive(Subcommand, Debug)]
pub enum MapCmd {
    /// Time-averaged number density [atoms/Å³] per voxel
    Density(GridCmd),
    /// Time-averaged speed |v| per voxel (needs frame velocities)
    Velocity(GridCmd),
    /// Time-averaged force magnitude |f| per voxel (needs frame forces)
    Force(GridCmd),
    /// Hard-sphere spatial occupancy map
    Radius(RadiusCmd),
    /// Cluster spatial distribution function (Qn-type, Kabsch alignment)
    Sdf(SdfCmd),
    /// Charge-density cluster SDF from QE pp.x cube files
    ChgSdf(ChgSdfCmd),
}

/// Voxel grid shared by density / velocity / force / radius.
#[derive(Args, Debug)]
pub struct GridKnobs {
    /// Grid points along the a axis
    #[arg(long, default_value = "50")]
    pub nx: usize,
    /// Grid points along the b axis
    #[arg(long, default_value = "50")]
    pub ny: usize,
    /// Grid points along the c axis
    #[arg(long, default_value = "50")]
    pub nz: usize,
    /// Only include these elements, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
}

#[derive(Args, Debug)]
pub struct GridCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub grid: GridKnobs,
}

#[derive(Args, Debug)]
pub struct RadiusCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub grid: GridKnobs,
    /// Hard-sphere radius cutoff [Å]
    #[arg(long, default_value = "0.7")]
    pub radius: f64,
}

/// Cluster identification shared by `sdf` and `chg-sdf`.
#[derive(Args, Debug)]
pub struct ClusterKnobs {
    /// Target Qn cluster level (0/1/2/3)
    #[arg(long, default_value = "3")]
    pub qn: u8,
    /// Network-former element
    #[arg(long, default_value = "P")]
    pub former: String,
    /// Ligand (bridging) element
    #[arg(long, default_value = "O")]
    pub ligand: String,
    /// Former-ligand bond cutoff [Å]
    #[arg(long, default_value = "2.4")]
    pub cutoff_fl: f64,
    /// Modifier element (e.g. Zn); omit if none
    #[arg(long)]
    pub modifier: Option<String>,
    /// Modifier-ligand cutoff [Å]
    #[arg(long, default_value = "2.8")]
    pub cutoff_ml: f64,
    /// RMSD warning threshold [Å]
    #[arg(long, default_value = "0.5")]
    pub rmsd_warn: f64,
}

#[derive(Args, Debug)]
pub struct SdfCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub cluster: ClusterKnobs,
    /// Voxel size [Å]
    #[arg(long, default_value = "0.1")]
    pub grid_res: f64,
    /// Gaussian broadening sigma [voxels]
    #[arg(long, default_value = "1.5")]
    pub sigma: f64,
    /// Grid boundary padding [Å]
    #[arg(long, default_value = "3.0")]
    pub padding: f64,
}

#[derive(Args, Debug)]
pub struct ChgSdfCmd {
    #[command(flatten)]
    pub cluster: ClusterKnobs,
    /// QE pp.x cube files (one per MD frame)
    #[arg(long = "cubes", num_args = 1..)]
    pub cube_files: Vec<PathBuf>,
    /// Output file stem
    #[arg(short, long)]
    pub output: Option<String>,
    /// Directory to write the cubes into (created if missing; default: current dir)
    #[arg(long, value_name = "DIR")]
    pub outdir: Option<PathBuf>,
    /// Sub-grid boundary margin [Å]
    #[arg(long = "chg-padding", default_value = "6.0")]
    pub chg_padding: f64,
    /// Parallel threads (default: all CPU cores)
    #[arg(long)]
    pub ncore: Option<usize>,
}

// ─── 分派 ────────────────────────────────────────────────────────────────────

pub fn wants_help(cmd: &MapCmd) -> bool {
    match cmd {
        MapCmd::Density(c) | MapCmd::Velocity(c) | MapCmd::Force(c) => c.common.input.is_empty(),
        MapCmd::Radius(c) => c.common.input.is_empty(),
        MapCmd::Sdf(c) => c.common.input.is_empty(),
        MapCmd::ChgSdf(c) => c.cube_files.is_empty(),
    }
}

pub fn print_help(cmd: &MapCmd) {
    match cmd {
        MapCmd::Density(_)  => help::print_cube_density(),
        MapCmd::Velocity(_) => help::print_cube_velocity(),
        MapCmd::Force(_)    => help::print_cube_force(),
        MapCmd::Radius(_)   => help::print_cube_radius(),
        MapCmd::Sdf(_)      => help::print_cube_sdf(),
        MapCmd::ChgSdf(_)   => help::print_cube_chg_sdf(),
    }
}

pub fn run(cmd: &MapCmd) -> Result<usize> {
    match cmd {
        MapCmd::Density(c)  => run_grid(c, CubeCliMode::Density),
        MapCmd::Velocity(c) => run_grid(c, CubeCliMode::Velocity),
        MapCmd::Force(c)    => run_grid(c, CubeCliMode::Force),
        MapCmd::Radius(c)   => run_radius(c),
        MapCmd::Sdf(c)      => run_sdf(c),
        MapCmd::ChgSdf(c)   => run_chg_sdf(c).map(|_| 0),
    }
}

/// Output stem: `<-o 或模式名>[_<输入 stem>]`.
///
/// 单输入时保持简洁名字；多输入时必须掺入 stem，否则第二个文件会盖掉第一个。
fn stem_for(output: Option<&str>, path: &Path, multi: bool) -> String {
    let base = output.unwrap_or("");
    let base = if base.is_empty() { String::new() } else { format!("{base}_") };
    if multi {
        format!("{base}{}", batch::label_of(path))
    } else {
        base.trim_end_matches('_').to_string()
    }
}

fn mode_name(m: &CubeCliMode) -> &'static str {
    match m {
        CubeCliMode::Density  => "density",
        CubeCliMode::Velocity => "velocity",
        CubeCliMode::Force    => "force",
        CubeCliMode::Radius   => "radius",
        CubeCliMode::Sdf      => "sdf",
        CubeCliMode::ChgSdf   => "chg_sdf",
    }
}

/// Shared driver: traverse inputs, skip failures, one grid file each.
///
/// `map` has no type selection, so its `Output` carries no label — only `--outdir` and
/// the `-o` stem, the latter already folded into `stem_for`.
fn drive(
    common: &CommonArgs,
    body: impl Fn(&Trajectory, &str, &batch::Output) -> Result<()>,
) -> Result<usize> {
    let out = common.out(None);
    out.prepare()?;

    let inputs = batch::expand_inputs(&common.input)?;
    common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let multi = inputs.len() > 1;
    let (_ok, failures) = batch::map_inputs(&inputs, |p| {
        let traj = common.load(p)?;
        let stem = stem_for(common.suffix(), p, multi);
        body(&traj, &stem, &out)
    });
    Ok(failures.len())
}

fn run_grid(c: &GridCmd, mode: CubeCliMode) -> Result<usize> {
    drive(&c.common, |traj, stem, out| {
        let params = CubeDensityParams {
            nx: c.grid.nx,
            ny: c.grid.ny,
            nz: c.grid.nz,
            elements: c.grid.elements.clone(),
            mode: mode.clone().into(),
        };
        let result = calc_cube_density(traj, &params)
            .ok_or_else(|| anyhow!("Cube calc failed (missing cell, velocities, or forces?)"))?;

        let name = mode_name(&mode);
        let file = if stem.is_empty() { format!("{name}.cube") } else { format!("{name}_{stem}.cube") };
        let path = out.join_str(&file);
        write_cube(&path, &result.cube)?;
        println!(
            "Cube ({name}) -> {path}  [{} frames, {} atoms]",
            result.n_frames, result.n_atoms
        );
        Ok(())
    })
}

fn run_radius(c: &RadiusCmd) -> Result<usize> {
    drive(&c.common, |traj, stem, out| {
        let params = CubeRadiusParams {
            nx: c.grid.nx,
            ny: c.grid.ny,
            nz: c.grid.nz,
            radius: c.radius,
            elements: c.grid.elements.clone(),
        };
        let result = calc_cube_radius(traj, &params)
            .ok_or_else(|| anyhow!("Cube radius calc failed (missing cell?)"))?;

        let file = if stem.is_empty() { "radius.cube".to_string() } else { format!("radius_{stem}.cube") };
        let path = out.join_str(&file);
        write_cube(&path, &result.cube)?;
        println!(
            "Cube (radius={:.3}Å) -> {path}  [{} frames, {} atoms]",
            params.radius, result.n_frames, result.n_atoms
        );
        Ok(())
    })
}

fn run_sdf(c: &SdfCmd) -> Result<usize> {
    drive(&c.common, |traj, stem, out| {
        let params = ClusterSdfParams {
            former: c.cluster.former.clone(),
            ligand: c.cluster.ligand.clone(),
            target_qn: c.cluster.qn,
            former_ligand_cutoff: c.cluster.cutoff_fl,
            modifier: c.cluster.modifier.clone(),
            modifier_cutoff: c.cluster.cutoff_ml,
            grid_res: c.grid_res,
            sigma: c.sigma,
            padding: c.padding,
            rmsd_warn_threshold: c.cluster.rmsd_warn,
        };
        let result = calc_cluster_sdf(traj, &params)
            .ok_or_else(|| anyhow!("No Q{} clusters found in trajectory", c.cluster.qn))?;

        let stem = if stem.is_empty() { "sdf" } else { stem };
        let multi_family = result.families.len() > 1;

        // 按签名排序,保证文件命名可复现
        let mut families: Vec<_> = result.families.values().collect();
        families.sort_by(|a, b| a.signature.cmp(&b.signature));

        let mut total_files = 0usize;
        for (fam_idx, family) in families.iter().enumerate() {
            let fam_prefix = if multi_family {
                format!("{stem}_fam{fam_idx}")
            } else {
                stem.to_string()
            };
            let mut labels: Vec<_> = family.grids.keys().collect();
            labels.sort();
            for label in &labels {
                write_cube(&out.join_str(&format!("{fam_prefix}_{label}.cube")), &family.grids[*label])?;
                total_files += 1;
            }
            println!(
                "Family {:?}  ({} clusters, RMSD mean={:.3} max={:.3} Å, {} warnings)  → {}_*.cube",
                family.signature,
                family.n_clusters,
                family.rmsd_stats.mean,
                family.rmsd_stats.max,
                family.rmsd_stats.n_warned,
                out.join_str(&fam_prefix),
            );
        }
        println!(
            "SDF Q{} done: {} frames, {} clusters total, {total_files} cube files written",
            c.cluster.qn, result.n_frames, result.n_clusters_total,
        );
        Ok(())
    })
}

/// `chg-sdf` takes `--cubes` instead of `-i`: it aggregates many cube files into ONE
/// SDF, the opposite of the batch model. Splitting it into a per-file command plus a
/// weighted cross-file average needs a new intermediate product format and is tracked
/// separately in `dev/plan.md`.
fn run_chg_sdf(c: &ChgSdfCmd) -> Result<()> {
    // 这个命令不走 CommonArgs(它吃 --cubes 而不是 -i),所以 Output 自己拼
    let out = batch::Output { dir: c.outdir.clone(), ..Default::default() };
    out.prepare()?;

    if let Some(n) = c.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let mut pairs = Vec::with_capacity(c.cube_files.len());
    for path in &c.cube_files {
        let (frame, chg) = read_cube_as_chg(path.to_str().unwrap_or_default())
            .map_err(|e| anyhow!("读取 {} 失败: {e}", path.display()))?;
        pairs.push((frame, chg));
    }

    let params = ChgSdfParams {
        former: c.cluster.former.clone(),
        ligand: c.cluster.ligand.clone(),
        target_qn: c.cluster.qn,
        former_ligand_cutoff: c.cluster.cutoff_fl,
        modifier: c.cluster.modifier.clone(),
        modifier_cutoff: c.cluster.cutoff_ml,
        padding: c.chg_padding,
        rmsd_warn_threshold: c.cluster.rmsd_warn,
    };
    let result = calc_chg_sdf(&pairs, &params)
        .ok_or_else(|| anyhow!("未找到 Q{} 团簇，请检查参数", c.cluster.qn))?;

    let stem = c.output.as_deref().unwrap_or("chg_sdf");
    let multi_family = result.families.len() > 1;
    let mut families: Vec<_> = result.families.values().collect();
    families.sort_by(|a, b| a.signature.cmp(&b.signature));

    let mut total_files = 0usize;
    for (fam_idx, family) in families.iter().enumerate() {
        let fam_stem = if multi_family { format!("{stem}_fam{fam_idx}") } else { stem.to_string() };
        let path = out.join_str(&format!("{}_Q{}.cube", fam_stem, c.cluster.qn));
        write_cube(&path, &family.cube)?;
        total_files += 1;
        println!(
            "Family {:?}  ({} clusters, RMSD mean={:.3} max={:.3} Å, {} warnings)  → {path}",
            family.signature,
            family.n_clusters,
            family.rmsd_stats.mean,
            family.rmsd_stats.max,
            family.rmsd_stats.n_warned,
        );
    }
    println!(
        "ChgSDF Q{} done: {} frames, {} clusters total, {total_files} cube files written",
        c.cluster.qn, result.n_frames, result.n_clusters_total,
    );
    Ok(())
}
