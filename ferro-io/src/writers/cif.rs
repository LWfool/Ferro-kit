use std::path::Path;
use std::collections::HashMap;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// Writes one structure as a CIF `data_` block in fractional coordinates
/// (`_atom_site_fract_*`) with P1 symmetry. Several frames are refused: CIF is written
/// one structure per file, since VESTA / Materials Studio show only a first block.
///
/// A cell-less frame is an error: CIF describes a crystal, and a non-periodic structure
/// belongs in xyz / extxyz (as with the POSCAR writer). It used to be written as a
/// Cartn-only block without cell parameters, which ferro itself could not read back.
pub fn write_cif(trajectory: &Trajectory, path: &Path) -> Result<()> {
    let path_ = path.display();
    // Checked before the file is created, so a refusal leaves nothing half-written
    let frame = super::single_frame(trajectory, "CIF")?;
    let cell = frame.cell.as_ref().context(
        "the frame has no cell; CIF requires a periodic frame (write a non-periodic \
         structure as .xyz or .extxyz)")?;
    let file = File::create(path).with_context(|| format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    // data block name
    let block_name = trajectory.metadata.source
        .as_deref()
        .map(sanitize_block_name)
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "structure".to_string());

    writeln!(w, "data_{block_name}")?;
    writeln!(w)?;

    let [a, b, c] = cell.lengths();
    let [alpha, beta, gamma] = cell.angles();

    writeln!(w, "_cell_length_a         {a:.6}")?;
    writeln!(w, "_cell_length_b         {b:.6}")?;
    writeln!(w, "_cell_length_c         {c:.6}")?;
    writeln!(w, "_cell_angle_alpha       {alpha:.4}")?;
    writeln!(w, "_cell_angle_beta        {beta:.4}")?;
    writeln!(w, "_cell_angle_gamma       {gamma:.4}")?;
    writeln!(w)?;
    writeln!(w, "_space_group_name_H-M   'P 1'")?;
    writeln!(w, "_space_group_IT_number  1")?;
    writeln!(w)?;

    writeln!(w, "loop_")?;
    writeln!(w, "_atom_site_label")?;
    writeln!(w, "_atom_site_type_symbol")?;
    writeln!(w, "_atom_site_fract_x")?;
    writeln!(w, "_atom_site_fract_y")?;
    writeln!(w, "_atom_site_fract_z")?;

    let mut counters: HashMap<String, usize> = HashMap::new();
    for atom in frame.atoms() {
        let label = if let Some(l) = atom.label {
            l.to_string()
        } else {
            let n = counters.entry(atom.element.to_string()).or_insert(0);
            *n += 1;
            format!("{}{n}", atom.element)
        };
        let frac = cell.cartesian_to_fractional(atom.position)?;
        writeln!(
            w,
            "{label:<8} {:<4} {:>10.6} {:>10.6} {:>10.6}",
            atom.element, frac.x, frac.y, frac.z
        )?;
    }

    writeln!(w)?;
    w.flush()?;
    Ok(())
}

fn sanitize_block_name(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::cif::read_cif;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn make_bcc_frame() -> Frame {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        frame
    }

    #[test]
    fn test_roundtrip_periodic() {
        let path = std::env::temp_dir().join("bcc_fe_rt.cif");
        let path_str = &path;

        let mut traj = Trajectory::from_frame(make_bcc_frame());
        traj.metadata.source = Some("BCC_Fe".to_string());
        write_cif(&traj, path_str).unwrap();

        let loaded = read_cif(path_str).unwrap();
        assert_eq!(loaded.n_frames(), 1);
        let frame = loaded.first().unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");

        let [a, ..] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }

    #[test]
    fn test_cell_less_frame_is_an_error() {
        let path = std::env::temp_dir().join("water_nocell.cif");
        let _ = std::fs::remove_file(&path);

        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.119)));
        let e = write_cif(&Trajectory::from_frame(frame), &path).unwrap_err();
        assert!(format!("{e:#}").contains("has no cell"), "{e:#}");
        assert!(format!("{e:#}").contains(".xyz"), "报错应指向 xyz：{e:#}");
        assert!(!path.exists(), "报错时不应留下半截文件");
    }
}
