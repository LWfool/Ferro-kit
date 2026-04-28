use std::path::Path;
use anyhow::{bail, Result};
use ferro_core::Trajectory;
use ferro_io::{self, LammpsUnits, *};

pub fn read_trajectory(path: &Path, lammps_units: LammpsUnits) -> Result<Trajectory> {
    let s = path.to_str().unwrap_or_default();
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();

    if upper.starts_with("POSCAR") {
        return Ok(read_poscar(s)?);
    }
    if upper.starts_with("CONTCAR") {
        return Ok(read_contcar(s)?);
    }

    match path.extension().and_then(|e| e.to_str()) {
        Some("xyz")                      => Ok(read_xyz(s)?),
        Some("pdb")                      => Ok(read_pdb(s)?),
        Some("cif")                      => Ok(read_cif(s)?),
        Some("extxyz")                   => Ok(read_extxyz(s)?),
        Some("lammps") | Some("data") | Some("lmp") => Ok(read_lammps_data(s)?),
        Some("dump") | Some("lammpstrj")             => Ok(read_lammps_dump(s, lammps_units)?),
        Some("inp")                      => Ok(read_cp2k_inp(s)?),
        Some("restart")                  => Ok(read_cp2k_restart(s)?),
        Some("in") | Some("qe")          => Ok(read_qe_input(s)?),
        Some(ext) => bail!("Unsupported input format: .{ext}"),
        None      => bail!("Cannot determine format (no extension): {s}"),
    }
}

pub fn write_trajectory(traj: &Trajectory, path: &Path, lammps_units: LammpsUnits) -> Result<()> {
    let s = path.to_str().unwrap_or_default();
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();

    if upper.starts_with("POSCAR") || upper.starts_with("CONTCAR") {
        return Ok(write_poscar(traj, s)?);
    }

    match path.extension().and_then(|e| e.to_str()) {
        Some("xyz")                      => Ok(write_xyz(traj, s)?),
        Some("pdb")                      => Ok(write_pdb(traj, s)?),
        Some("cif")                      => Ok(write_cif(traj, s)?),
        Some("extxyz")                   => Ok(write_extxyz(traj, s)?),
        Some("lammps") | Some("data") | Some("lmp") => Ok(write_lammps_data(traj, s)?),
        Some("dump") | Some("lammpstrj") => Ok(write_lammps_dump(traj, s, lammps_units)?),
        Some("in") | Some("qe")          => Ok(write_qe_input(traj, s)?),
        Some(ext) => bail!("Unsupported output format: .{ext}"),
        None      => bail!("Cannot determine format (no extension): {s}"),
    }
}

/// Returns a human-readable list of supported read/write formats.
pub fn supported_formats() -> &'static str {
    "Read:  xyz, pdb, cif, extxyz, lammps/data, dump/lammpstrj, inp, restart, in/qe, POSCAR, CONTCAR
Write: xyz, pdb, cif, extxyz, lammps/data, dump/lammpstrj, in/qe, POSCAR"
}
