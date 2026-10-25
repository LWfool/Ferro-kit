use std::path::Path;
use std::collections::{HashMap, HashSet};
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
    let block_name = trajectory.metadata.title
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

    // `_atom_site_label` is unique (review R2-7): an automatic label skips every name an
    // explicit label holds, anywhere in the frame, and a repeated explicit label gets
    // `_2`, `_3`, ...
    let explicit: HashSet<&str> = frame.atoms().filter_map(|a| a.label).collect();
    let mut used: HashSet<String> = HashSet::new();
    let mut counters: HashMap<String, usize> = HashMap::new();
    for atom in frame.atoms() {
        let label = if let Some(l) = atom.label {
            (1..).map(|k| if k == 1 { l.to_string() } else { format!("{l}_{k}") })
                .find(|c| !used.contains(c) && (c == l || !explicit.contains(c.as_str())))
                .expect("an unbounded range finds a free name")
        } else {
            let n = counters.entry(atom.element.to_string()).or_insert(0);
            loop {
                *n += 1;
                let c = format!("{}{n}", atom.element);
                if !used.contains(&c) && !explicit.contains(c.as_str()) { break c; }
            }
        };
        used.insert(label.clone());
        let frac = cell.cartesian_to_fractional(atom.position)?;
        writeln!(
            w,
            "{:<8} {:<4} {:>10.6} {:>10.6} {:>10.6}",
            quote(&label), atom.element, frac.x, frac.y, frac.z
        )?;
    }

    writeln!(w)?;
    w.flush()?;
    Ok(())
}

/// A CIF value as one token: quoted when it holds a blank or starts with a character
/// that would make it something else (a tag, a comment, a quoted or text field)
fn quote(s: &str) -> String {
    let special = s.starts_with(['_', '#', '$', '\'', '"', ';', '[', ']']);
    if !s.contains(char::is_whitespace) && !special {
        s.to_string()
    } else if s.contains('\'') {
        format!("\"{s}\"")
    } else {
        format!("'{s}'")
    }
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
    use ferro_core::{Atom, Frame, Trajectory};
    use nalgebra::Vector3;

    /// `_atom_site_label` must be unique and one token. An automatic label used to take a
    /// name an explicit label held, a repeated explicit label was written twice, and a
    /// label with a space shifted the loop's columns.
    #[test]
    fn test_labels_are_unique_and_quoted() {
        let cell = ferro_core::Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        let mut add = |el: &str, label: Option<&str>, x: f64| {
            let mut a = Atom::new(el, Vector3::new(x, 0.0, 0.0));
            a.label = label.map(Into::into);
            frame.add_atom(a);
        };
        add("Fe", None, 0.0);
        add("Fe", Some("Fe1"), 1.0);
        add("O", Some("O a"), 2.0);
        add("O", Some("O a"), 3.0);
        let path = std::env::temp_dir().join("ferro_cif_labels.cif");
        write_cif(&Trajectory::from_frame(frame), &path).unwrap();
        let back = read_cif(&path).unwrap();
        let labels: Vec<_> = back.frames[0].atoms().map(|a| a.label.map(str::to_string)).collect();
        assert_eq!(labels, [Some("Fe2"), Some("Fe1"), Some("O a"), Some("O a_2")].map(|l| l.map(Into::into)));
    }

    #[test]
    fn test_roundtrip_periodic() {
        let path = std::env::temp_dir().join("bcc_fe_rt.cif");
        let path_str = &path;

        let mut traj = Trajectory::from_frame(crate::testutil::bcc_fe());
        traj.metadata.title = Some("BCC_Fe".to_string());
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
