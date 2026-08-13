//! Argument groups shared by several subcommands, flattened in with
//! `#[command(flatten)]` so each subcommand declares them once instead of copying
//! the fields.

use anyhow::{anyhow, Result};
use clap::Args;
use ferro_analysis::GroupBy;
use ferro_core::Trajectory;
use ferro_io::LammpsUnits;
use std::path::{Path, PathBuf};

use crate::io_dispatch::read_trajectory_tail;

/// Input/output and run-wide options carried by every analysis subcommand.
#[derive(Args, Clone, Debug)]
pub struct CommonArgs {
    /// Input file(s); glob patterns allowed — quote them so the shell leaves them alone
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    pub input: Vec<PathBuf>,

    /// Output name suffix: results go to <command>[_<table>][_<label>]_<suffix>.csv
    #[arg(short, long, value_name = "SUFFIX")]
    pub output: Option<String>,

    /// Directory to write every product into (created if missing; default: current dir)
    #[arg(long, value_name = "DIR")]
    pub outdir: Option<PathBuf>,

    /// Use only the last N frames (skip equilibration)
    #[arg(long)]
    pub last_n: Option<usize>,

    /// Parallel threads (default: all CPU cores)
    #[arg(long)]
    pub ncore: Option<usize>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    pub metal_units: bool,
}

impl CommonArgs {
    pub fn units(&self) -> LammpsUnits {
        if self.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real }
    }

    /// Applies `--ncore` to the global rayon pool. Call once, before any analysis.
    pub fn init_threads(&self) {
        if let Some(n) = self.ncore {
            rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
        }
    }

    /// Reads one trajectory and applies `--last-n`.
    pub fn load(&self, path: &Path) -> Result<Trajectory> {
        read_trajectory_tail(path, self.units(), self.last_n)
    }

    pub fn suffix(&self) -> Option<&str> {
        self.output.as_deref()
    }

    /// Builds the naming/placement bundle for this run.
    ///
    /// `label` says what was analysed and lands in the file name; pass `None` for the
    /// modes that have no type selection (`sq`, `net`, `map`).
    pub fn out(&self, label: Option<String>) -> crate::batch::Output {
        crate::batch::Output {
            dir: self.outdir.clone(),
            label,
            suffix: self.output.clone(),
        }
    }
}

/// Type selection for the pair/triplet analyses.
///
/// Flattened into `gr` and `angle` only — **not `sq`**, whose partials are always all
/// written.
///
/// Two mutually exclusive groups: `-a/-b/-c` pick chemical elements, `-x/-y/-z` pick
/// site labels. Position carries meaning — the first slot is the centre for a pair and
/// end A for a triplet — so the order the caller wrote is preserved downstream rather
/// than being re-sorted by atomic number.
#[derive(Args, Clone, Debug, Default)]
pub struct SelectArgs {
    /// Centre element (gr) or end atom A (angle); requires -b
    #[arg(short = 'a', long)]
    pub atom_a: Option<String>,

    /// Neighbour element (gr) or centre atom B (angle); requires -a
    #[arg(short = 'b', long)]
    pub atom_b: Option<String>,

    /// End atom C by element (angle only); requires -a -b
    #[arg(short = 'c', long)]
    pub atom_c: Option<String>,

    /// Centre site label (gr) or end atom A (angle); requires -y
    #[arg(short = 'x', long)]
    pub label_x: Option<String>,

    /// Neighbour site label (gr) or centre atom B (angle); requires -x
    #[arg(short = 'y', long)]
    pub label_y: Option<String>,

    /// End atom C by site label (angle only); requires -x -y
    #[arg(short = 'z', long)]
    pub label_z: Option<String>,
}

impl SelectArgs {
    /// Resolves the grouping mode and the three slots, rejecting a mixed selection.
    pub fn resolve(&self) -> Result<(GroupBy, [Option<String>; 3])> {
        let by_elem = [&self.atom_a, &self.atom_b, &self.atom_c];
        let by_label = [&self.label_x, &self.label_y, &self.label_z];
        let n_elem = by_elem.iter().filter(|x| x.is_some()).count();
        let n_label = by_label.iter().filter(|x| x.is_some()).count();

        if n_elem > 0 && n_label > 0 {
            return Err(anyhow!(
                "-a/-b/-c (by element) and -x/-y/-z (by label) are mutually exclusive; use one group"
            ));
        }
        let slots = if n_label > 0 { by_label } else { by_elem };
        let mode = if n_label > 0 { GroupBy::Label } else { GroupBy::Element };
        Ok((mode, [slots[0].clone(), slots[1].clone(), slots[2].clone()]))
    }

    /// Pair selection for gr: either both members or neither.
    ///
    /// `sq` used to share this and no longer does — it has no type selection at all.
    pub fn resolve_pair(&self) -> Result<(GroupBy, Option<(String, String)>)> {
        let (mode, slots) = self.resolve()?;
        if slots[2].is_some() {
            return Err(anyhow!("-c / -z applies to the angle command only"));
        }
        match (&slots[0], &slots[1]) {
            (Some(a), Some(b)) => Ok((mode, Some((a.clone(), b.clone())))),
            (None, None) => Ok((mode, None)),
            _ => Err(anyhow!(
                "pair filter needs both members: -a CENTRE -b NEIGHBOUR (by element) \
                 or -x CENTRE -y NEIGHBOUR (by label)"
            )),
        }
    }

    /// Triplet selection for angle: all three or none.
    pub fn resolve_triplet(&self) -> Result<(GroupBy, [Option<String>; 3])> {
        let (mode, slots) = self.resolve()?;
        let n_given = slots.iter().filter(|x| x.is_some()).count();
        if n_given > 0 && n_given < 3 {
            return Err(anyhow!(
                "angle filter needs all three: -a END_A -b CENTRE -c END_C (by element) \
                 or -x END_A -y CENTRE -z END_C (by label)"
            ));
        }
        Ok((mode, slots))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sel(a: Option<&str>, b: Option<&str>, c: Option<&str>, x: Option<&str>) -> SelectArgs {
        SelectArgs {
            atom_a: a.map(str::to_string),
            atom_b: b.map(str::to_string),
            atom_c: c.map(str::to_string),
            label_x: x.map(str::to_string),
            ..Default::default()
        }
    }

    #[test]
    fn test_mixed_groups_rejected() {
        let s = sel(Some("P"), None, None, Some("P_0"));
        assert!(s.resolve().is_err(), "元素组与标签组不能混用");
    }

    #[test]
    fn test_pair_needs_both_members() {
        assert!(sel(Some("P"), None, None, None).resolve_pair().is_err());
        assert!(sel(None, None, None, None).resolve_pair().is_ok(), "两个都不给是合法的");
        let (_, pair) = sel(Some("P"), Some("O"), None, None).resolve_pair().unwrap();
        assert_eq!(pair, Some(("P".into(), "O".into())), "顺序原样保留,中心在前");
    }

    #[test]
    fn test_pair_rejects_third_slot() {
        assert!(sel(Some("O"), Some("P"), Some("O"), None).resolve_pair().is_err());
    }

    #[test]
    fn test_triplet_is_all_or_nothing() {
        assert!(sel(Some("O"), Some("P"), None, None).resolve_triplet().is_err());
        assert!(sel(Some("O"), Some("P"), Some("O"), None).resolve_triplet().is_ok());
        assert!(sel(None, None, None, None).resolve_triplet().is_ok());
    }

    #[test]
    fn test_label_group_switches_mode() {
        let s = SelectArgs { label_x: Some("P_0".into()), label_y: Some("O_f".into()), ..Default::default() };
        let (mode, _) = s.resolve_pair().unwrap();
        assert!(matches!(mode, GroupBy::Label));
    }
}
