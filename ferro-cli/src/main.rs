fn main() {
    println!(
        r#"ferro — Computational Chemistry Toolkit  v{}

Available commands:
  fe-convert   Format conversion  (xyz, pdb, cif, poscar, extxyz, lammps, qe, ...)
  fe-info      Display structure / trajectory information
  fe-job       Generate QC input files  (gaussian, ...)
  fe-traj      Structural analysis      -m gr | sq | msd | angle
  fe-corr      Correlation functions    -m vacf | rotcorr | vanhove
  fe-cube      Spatial density maps     -m density | velocity | force
  fe-network   Glass network analysis   CN, FO/NBO/BO/OBO, Qn  (--Former-Ligand=cutoff)

Usage:
  <command> --help           Show all flags
  <command> -m <mode>        Show mode-specific help (fe-traj / fe-corr / fe-cube)
  <command> -m <mode> -i ...  Run analysis"#,
        env!("CARGO_PKG_VERSION")
    );
}
