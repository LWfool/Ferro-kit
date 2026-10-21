// One public path per item (review A-12): everything is reached from the crate root
mod format;
mod notice;
mod readers;
mod writers;

pub use format::{
    is_deepmd_system, read, read_as, read_charge_grid, write, write_numbered, Format, ReadOptions, ReadReport,
    WriteOptions,
};

pub use notice::{Level, Notice};

#[cfg(test)]
pub(crate) mod testutil;

pub use readers::{
    read_aimd_as, read_aimd, read_cp2k_md,
    read_cp2k_sp,
    read_vasp_outcar,
    read_vasprun, sniff, AimdFormat, AimdStats,
    read_deepmd_npy,
    read_cif, read_pdb, read_xyz,
    read_poscar,
    read_chgcar,
    read_cube_as_chg,
    is_extxyz, read_extxyz,
    read_lammps_data, read_lammps_dump, apply_type_map,
    read_cp2k_inp, read_cp2k_restart,
    read_qe_out,
    read_qe_input, parse_qe_input_parts, write_namelists, Namelist, QeInputParts,
    AtomStyle, LammpsUnits,
};
pub use writers::{
    set_bounds, write_deepmd_npy, write_deepmd_npy_bounds, write_deepmd_npy_sets,
    write_cif, write_pdb, write_xyz,
    write_poscar, poscar_string,
    write_extxyz, write_extxyz_with, StressKey,
    write_lammps_data, write_lammps_dump, dump_name,
    qe_structure_cards, cp2k_structure_sections,
    write_cube,
    write_table, TableFormat,
    index_width, indexed_path,
};
