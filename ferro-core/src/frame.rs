//! Single-frame structure (analogous to ASE's `Atoms` object).

use std::sync::Arc;

use nalgebra::{Matrix3, Vector3};

use crate::atom::{Atom, AtomRef, Site};
use crate::cell::Cell;
use crate::error::{ChemError, Result};

/// One structure: the core unit of ferro.
///
/// - Non-periodic (molecule): `cell = None`, `pbc = [false; 3]`
/// - 3-D periodic (crystal): `cell = Some(…)`, `pbc = [true; 3]`
/// - 2-D periodic (surface): `cell = Some(…)`, `pbc = [true, true, false]`
///
/// An NPT trajectory has a different box in every frame; each frame of a
/// [`crate::trajectory::Trajectory`] holding its own `cell` handles that without
/// anything special.
///
/// **Storage** (review B15). What an atom is — element, label, mass, magnetic
/// moment — sits in a table of [`Site`]s behind an `Arc`, which consecutive frames
/// of one system share ([`crate::Trajectory::add_frame`] links them); what changes
/// from frame to frame — positions, charges, forces, velocities — sits in arrays of
/// the frame's own. A frame costs 24 bytes per atom for its positions instead of a
/// whole `Atom` with two strings. The per-atom arrays are private so their lengths
/// cannot drift from the atom count (review B5): [`Frame::set_forces`] and
/// [`Frame::set_velocities`] refuse a wrong length.
#[derive(Debug, Clone)]
pub struct Frame {
    sites: Arc<Vec<Site>>,
    positions: Vec<Vector3<f64>>,
    /// Per-atom charge \[e\]; `None` when no atom has one
    charges: Option<Vec<Option<f64>>>,
    /// Per-atom force \[eV/Å\], in atom order
    forces: Option<Vec<Vector3<f64>>>,
    /// Per-atom velocity \[Å/fs\], in atom order
    velocities: Option<Vec<Vector3<f64>>>,
    /// Lattice; `None` for a non-periodic system
    pub cell: Option<Cell>,
    /// Whether each direction is periodic
    pub pbc: [bool; 3],

    // ── Results ──────────────────────────────────────────────────────────────
    /// Total energy \[eV\]
    pub energy: Option<f64>,
    /// Stress tensor (eV/Å³), row-major, **positive = compression**.
    ///
    /// This is the sign CP2K, VASP and QE print, and the opposite of ASE's
    /// (positive = tension); readers of ASE-flavoured formats negate on the way
    /// in. The virial follows from it without a sign change: `virial = stress *
    /// V` (eV), which is what DeePMD, QUIP and GPUMD call `virial`.
    ///
    /// Voigt order is left to the caller — nothing in this crate flattens it.
    pub stress: Option<Matrix3<f64>>,
    /// Instantaneous ionic temperature (K), as an MD engine reports it per step.
    ///
    /// `None` for everything that is not an MD frame — a relaxation, a POSCAR, a
    /// CIF. It lives here rather than in a parallel array because it is a
    /// per-frame quantity: `select`/`stride` must carry it along, and a side
    /// channel would silently shift out of step with the frames.
    ///
    /// CP2K and VASP OUTCAR print it directly; vasprun.xml does not, so that
    /// reader derives it from the ionic kinetic energy and says so.
    pub temperature: Option<f64>,
    /// The MD step number the engine gave this frame, when it prints one.
    ///
    /// Kept because `collect` stitches restart segments without de-duplicating
    /// them: overlapping frames are genuinely identical, and the step number is
    /// what makes that overlap visible instead of silently doubling a run.
    /// CP2K and VASP OUTCAR number their ionic steps; vasprun.xml does not.
    pub step: Option<i64>,
}

impl Frame {
    /// An empty, non-periodic frame.
    pub fn new() -> Self {
        Self {
            sites: Arc::new(Vec::new()),
            positions: Vec::new(),
            charges: None,
            forces: None,
            velocities: None,
            cell: None,
            pbc: [false; 3],
            energy: None,
            stress: None,
            temperature: None,
            step: None,
        }
    }

    /// A periodic frame (the usual entry for crystals and surfaces).
    pub fn with_cell(cell: Cell, pbc: [bool; 3]) -> Self {
        Self {
            cell: Some(cell),
            pbc,
            ..Self::new()
        }
    }

    // ── Atoms ────────────────────────────────────────────────────────────────

    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Atom `index`, read out of the site table and the arrays.
    pub fn atom(&self, index: usize) -> AtomRef<'_> {
        let s = &self.sites[index];
        AtomRef {
            element: &s.element,
            position: self.positions[index],
            label: s.label.as_deref(),
            mass: s.mass,
            magmom: s.magmom,
            charge: self.charges.as_ref().and_then(|c| c[index]),
        }
    }

    /// Atom `index`, or `None` past the last one.
    pub fn get_atom(&self, index: usize) -> Option<AtomRef<'_>> {
        (index < self.n_atoms()).then(|| self.atom(index))
    }

    /// Every atom, in index order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomRef<'_>> + DoubleEndedIterator + Clone + '_ {
        (0..self.n_atoms()).map(move |i| self.atom(i))
    }

    /// Appends an atom.
    ///
    /// # Panics
    /// When the frame already carries forces or velocities: they would be one short.
    /// Add the atoms first, then set them.
    pub fn add_atom(&mut self, atom: Atom) {
        assert!(
            self.forces.is_none() && self.velocities.is_none(),
            "add_atom after set_forces / set_velocities would leave them one atom short"
        );
        let Atom { element, position, label, mass, magmom, charge } = atom;
        if charge.is_some() && self.charges.is_none() {
            self.charges = Some(vec![None; self.n_atoms()]);
        }
        if let Some(c) = &mut self.charges {
            c.push(charge);
        }
        Arc::make_mut(&mut self.sites).push(Site { element, label, mass, magmom });
        self.positions.push(position);
    }

    /// Replaces every atom; forces and velocities are dropped with the old atoms.
    pub fn set_atoms(&mut self, atoms: impl IntoIterator<Item = Atom>) {
        self.sites = Arc::new(Vec::new());
        self.positions.clear();
        self.charges = None;
        self.forces = None;
        self.velocities = None;
        for a in atoms {
            self.add_atom(a);
        }
    }

    /// The site table: element, label, mass and magnetic moment of every atom.
    pub fn sites(&self) -> &[Site] {
        &self.sites
    }

    /// The site table, for changing an element or a label. The table is copied
    /// first when other frames share it.
    pub fn sites_mut(&mut self) -> &mut [Site] {
        Arc::make_mut(&mut self.sites).as_mut_slice()
    }

    /// Whether two frames share one site table (so their atoms are the same by construction).
    pub fn shares_sites_with(&self, other: &Frame) -> bool {
        Arc::ptr_eq(&self.sites, &other.sites)
    }

    /// Makes this frame use `other`'s site table when the two are equal, so a
    /// trajectory stores it once.
    pub fn share_sites_with(&mut self, other: &Frame) {
        if !self.shares_sites_with(other) && self.sites == other.sites {
            self.sites = Arc::clone(&other.sites);
        }
    }

    /// Cartesian positions of all atoms, in index order.
    pub fn positions(&self) -> &[Vector3<f64>] {
        &self.positions
    }

    pub fn positions_mut(&mut self) -> &mut [Vector3<f64>] {
        &mut self.positions
    }

    /// Element symbol of atom `index`.
    pub fn element(&self, index: usize) -> &str {
        &self.sites[index].element
    }

    /// Sets the charge of atom `index`.
    pub fn set_charge(&mut self, index: usize, charge: Option<f64>) {
        let n = self.n_atoms();
        self.charges.get_or_insert_with(|| vec![None; n])[index] = charge;
    }

    /// Per-atom forces \[eV/Å\], in atom order.
    pub fn forces(&self) -> Option<&[Vector3<f64>]> {
        self.forces.as_deref()
    }

    pub fn forces_mut(&mut self) -> Option<&mut [Vector3<f64>]> {
        self.forces.as_deref_mut()
    }

    /// Sets the forces; one per atom, or an error.
    pub fn set_forces(&mut self, forces: Option<Vec<Vector3<f64>>>) -> Result<()> {
        self.check_len("forces", forces.as_ref())?;
        self.forces = forces;
        Ok(())
    }

    /// Per-atom velocities \[Å/fs\], in atom order.
    pub fn velocities(&self) -> Option<&[Vector3<f64>]> {
        self.velocities.as_deref()
    }

    pub fn velocities_mut(&mut self) -> Option<&mut [Vector3<f64>]> {
        self.velocities.as_deref_mut()
    }

    /// Sets the velocities; one per atom, or an error.
    pub fn set_velocities(&mut self, velocities: Option<Vec<Vector3<f64>>>) -> Result<()> {
        self.check_len("velocities", velocities.as_ref())?;
        self.velocities = velocities;
        Ok(())
    }

    fn check_len(&self, what: &str, v: Option<&Vec<Vector3<f64>>>) -> Result<()> {
        match v {
            Some(v) if v.len() != self.n_atoms() => Err(ChemError::ValidationError(format!(
                "{} {what} for {} atoms", v.len(), self.n_atoms()
            ))),
            _ => Ok(()),
        }
    }

    /// The frame with its atoms in the order `order` gives: new atom `k` is old atom
    /// `order[k]`. Every per-atom array moves with its atoms.
    ///
    /// # Panics
    /// When `order` is not a permutation of `0..n_atoms`.
    pub fn permuted(&self, order: &[usize]) -> Frame {
        let n = self.n_atoms();
        let mut seen = vec![false; n];
        assert!(
            order.len() == n && order.iter().all(|&i| i < n && !std::mem::replace(&mut seen[i], true)),
            "permuted: the order is not a permutation of 0..{n}"
        );
        let pick = |v: &[Vector3<f64>]| order.iter().map(|&i| v[i]).collect::<Vec<_>>();
        Frame {
            sites: Arc::new(order.iter().map(|&i| self.sites[i].clone()).collect()),
            positions: pick(&self.positions),
            charges: self.charges.as_ref().map(|c| order.iter().map(|&i| c[i]).collect()),
            forces: self.forces.as_deref().map(pick),
            velocities: self.velocities.as_deref().map(pick),
            ..self.clone_without_atoms()
        }
    }

    /// Keeps the first `n` atoms (and their forces, velocities and charges).
    pub fn truncate_atoms(&mut self, n: usize) {
        if n >= self.n_atoms() {
            return;
        }
        Arc::make_mut(&mut self.sites).truncate(n);
        self.positions.truncate(n);
        for v in [&mut self.forces, &mut self.velocities].into_iter().flatten() {
            v.truncate(n);
        }
        if let Some(c) = &mut self.charges {
            c.truncate(n);
        }
    }

    /// The frame-level fields, with no atoms.
    fn clone_without_atoms(&self) -> Frame {
        Frame {
            sites: Arc::new(Vec::new()),
            positions: Vec::new(),
            charges: None,
            forces: None,
            velocities: None,
            cell: self.cell.clone(),
            pbc: self.pbc,
            energy: self.energy,
            stress: self.stress,
            temperature: self.temperature,
            step: self.step,
        }
    }

    /// The virial (eV): `stress * V`, no sign change (see [`Frame::stress`]); `None`
    /// without a stress. A stress without a cell, or with a zero-volume one, has no
    /// virial and is an error. The inverse is [`Cell::stress_from_virial`].
    pub fn virial(&self) -> Result<Option<Matrix3<f64>>> {
        let Some(s) = &self.stress else { return Ok(None) };
        let cell = self.cell.as_ref().ok_or_else(|| ChemError::ValidationError(
            "a stress without a cell has no virial: the virial needs the volume".into(),
        ))?;
        Ok(Some(s * cell.nonzero_volume()?))
    }

    /// Element symbols of all atoms, in index order.
    pub fn symbols(&self) -> Vec<&str> {
        self.sites.iter().map(|s| s.element.as_str()).collect()
    }

    /// Distinct element symbols in order of first appearance (KIND / type ids / the POSCAR element line).
    pub fn unique_elements(&self) -> Vec<String> {
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for s in self.sites.iter() {
            if seen.insert(s.element.as_str()) {
                out.push(s.element.clone());
            }
        }
        out
    }

    /// Number of atoms of `element` in this frame.
    pub fn count_element(&self, element: &str) -> usize {
        self.sites.iter().filter(|s| s.element == element).count()
    }

    // ── Geometry ─────────────────────────────────────────────────────────────

    /// Total mass \[amu\].
    pub fn total_mass(&self) -> f64 {
        self.sites.iter().map(|s| s.effective_mass()).sum()
    }

    /// Mass-weighted centre \[Å\].
    pub fn center_of_mass(&self) -> Vector3<f64> {
        if self.positions.is_empty() {
            return Vector3::zeros();
        }
        let total = self.total_mass();
        self.sites
            .iter()
            .zip(&self.positions)
            .map(|(s, p)| p * s.effective_mass())
            .fold(Vector3::zeros(), |acc, v| acc + v)
            / total
    }

    /// Unweighted centre \[Å\].
    pub fn geometric_center(&self) -> Vector3<f64> {
        if self.positions.is_empty() {
            return Vector3::zeros();
        }
        let sum = self.positions.iter().fold(Vector3::zeros(), |acc, v| acc + v);
        sum / self.positions.len() as f64
    }

    // ── Transformations ──────────────────────────────────────────────────────

    /// Moves every atom by `disp`.
    pub fn translate(&mut self, disp: Vector3<f64>) {
        for p in &mut self.positions {
            *p += disp;
        }
    }

    /// Moves the centre of mass to the origin.
    pub fn center(&mut self) {
        let com = self.center_of_mass();
        self.translate(-com);
    }

    // ── Periodicity ──────────────────────────────────────────────────────────

    /// Whether any direction is periodic.
    pub fn is_periodic(&self) -> bool {
        self.pbc.iter().any(|&p| p)
    }

    /// Folds atoms outside the box back into `[0, 1)` along the periodic directions;
    /// atoms inside, and the non-periodic directions, are left alone.
    ///
    /// An atom outside moves by a whole number of lattice vectors, without a
    /// fractional -> Cartesian round trip: the round trip would leave a 1e-15
    /// rounding on every atom, inside ones included. "Outside" allows `WRAP_TOL`, so
    /// an atom at exactly 0 whose fraction comes out -1e-17 is not sent to the far
    /// side. Does nothing without a cell or with a singular one.
    pub fn wrap_all(&mut self) {
        const WRAP_TOL: f64 = 1e-10;
        let Some(cell) = self.cell.as_ref() else { return };
        let pbc = self.pbc;
        for p in &mut self.positions {
            let Ok(f) = cell.cartesian_to_fractional(*p) else { return };
            let shift = Vector3::from_fn(|i, _| {
                let out = f[i] < -WRAP_TOL || f[i] >= 1.0 + WRAP_TOL;
                if pbc[i] && out { -f[i].floor() } else { 0.0 }
            });
            if shift != Vector3::zeros() {
                *p += cell.fractional_to_cartesian(shift);
            }
        }
    }
}

impl Default for Frame {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;

    fn three_atoms() -> Frame {
        let mut f = Frame::new();
        for (e, x) in [("O", 0.0), ("H", 1.0), ("Fe", 2.0)] {
            let mut a = Atom::new(e, Vector3::new(x, 0.0, 0.0));
            a.label = Some(format!("{e}_x"));
            f.add_atom(a);
        }
        f
    }

    /// Per-atom arrays cannot drift from the atom count (review B5).
    #[test]
    fn forces_and_velocities_of_the_wrong_length_are_refused() {
        let mut f = three_atoms();
        let err = f.set_forces(Some(vec![Vector3::zeros(); 2])).unwrap_err().to_string();
        assert!(err.contains("2 forces for 3 atoms"), "{err}");
        assert!(f.forces().is_none());
        f.set_velocities(Some(vec![Vector3::zeros(); 3])).unwrap();
        assert_eq!(f.velocities().unwrap().len(), 3);
    }

    /// A permutation moves every per-atom quantity with its atom.
    #[test]
    fn permuted_moves_everything_with_its_atom() {
        let mut f = three_atoms();
        f.set_charge(1, Some(0.5));
        f.set_forces(Some((0..3).map(|i| Vector3::new(i as f64, 0.0, 0.0)).collect())).unwrap();
        let g = f.permuted(&[2, 0, 1]);
        assert_eq!(g.symbols(), ["Fe", "O", "H"]);
        assert_eq!(g.atom(2).label, Some("H_x"));
        assert_eq!(g.atom(2).charge, Some(0.5));
        assert_eq!(g.atom(0).charge, None);
        assert_eq!(g.forces().unwrap()[0].x, 2.0);
        assert_eq!(g.positions()[1].x, 0.0);
    }

    /// Stress and virial convert into each other, and a zero volume is refused both ways.
    #[test]
    fn virial_and_stress_round_trip_and_refuse_a_zero_volume() {
        let cell = Cell::from_lengths_angles(2.0, 3.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::with_cell(cell.clone(), [true; 3]);
        assert_eq!(f.virial().unwrap(), None);
        let s = Matrix3::new(1.0, 0.5, 0.0, 0.5, -2.0, 0.0, 0.0, 0.0, 3.0);
        f.stress = Some(s);
        let v = f.virial().unwrap().unwrap();
        assert_eq!(v, s * 24.0);
        assert_eq!(cell.stress_from_virial(&v).unwrap(), s);

        let flat = Cell::from_matrix(Matrix3::new(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0));
        assert!(flat.stress_from_virial(&v).unwrap_err().to_string().contains("volume is 0"));
        f.cell = Some(flat);
        assert!(f.virial().is_err());
        f.cell = None;
        assert!(f.virial().unwrap_err().to_string().contains("without a cell"));
    }

    /// Equal site tables are shared, and changing one frame's sites leaves the other's.
    #[test]
    fn site_tables_are_shared_and_copied_on_write() {
        let a = three_atoms();
        let mut b = three_atoms();
        assert!(!b.shares_sites_with(&a));
        b.share_sites_with(&a);
        assert!(b.shares_sites_with(&a));
        b.sites_mut()[0].element = "N".into();
        assert!(!b.shares_sites_with(&a));
        assert_eq!((a.element(0), b.element(0)), ("O", "N"));
    }

    #[test]
    fn test_new_frame() {
        let f = Frame::new();
        assert_eq!(f.n_atoms(), 0);
        assert!(!f.is_periodic());
    }

    #[test]
    fn test_add_and_count() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        f.add_atom(Atom::new("O",  Vector3::new(2.0, 0.0, 0.0)));
        f.add_atom(Atom::new("O",  Vector3::new(0.0, 2.0, 0.0)));
        assert_eq!(f.n_atoms(), 3);
        assert_eq!(f.count_element("O"), 2);
    }

    #[test]
    fn test_unique_elements_first_appearance() {
        let mut f = Frame::new();
        for e in ["O", "H", "O", "Fe", "H", "O"] {
            f.add_atom(Atom::new(e, Vector3::zeros()));
        }
        // 去重并保持首次出现顺序
        assert_eq!(f.unique_elements(), vec!["O", "H", "Fe"]);
        assert!(Frame::new().unique_elements().is_empty());
    }

    #[test]
    fn test_center_of_mass() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        f.add_atom(Atom::new("C", Vector3::new(2.0, 0.0, 0.0)));
        let com = f.center_of_mass();
        assert!((com.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_periodic_frame() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let f = Frame::with_cell(cell, [true; 3]);
        assert!(f.is_periodic());
        assert!(f.cell.is_some());
    }

    #[test]
    fn test_wrap_all() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.add_atom(Atom::new("Fe", Vector3::new(11.0, -1.0, 5.0)));
        f.wrap_all();
        let pos = f.atom(0).position;
        assert!(pos.x >= 0.0 && pos.x < 10.0);
        assert!(pos.y >= 0.0 && pos.y < 10.0);
    }

    #[test]
    fn test_wrap_all_respects_pbc_per_axis_and_leaves_inside_atoms_untouched() {
        // 三斜，z 非周期（slab）
        let cell = Cell::from_lengths_angles(9.0, 11.0, 8.0, 80.0, 105.0, 70.0).unwrap();
        let mut f = Frame::with_cell(cell.clone(), [true, true, false]);
        let inside = cell.fractional_to_cartesian(Vector3::new(0.3, 0.7, 0.2)) + Vector3::new(1e-13, 0.0, 0.0);
        let at_zero = Vector3::new(0.0, 0.0, 0.0);
        let out = cell.fractional_to_cartesian(Vector3::new(1.25, -0.4, 1.6));
        for p in [inside, at_zero, out] { f.add_atom(Atom::new("Si", p)); }
        f.wrap_all();
        // 盒内原子逐位不动（不做分数往返）
        assert_eq!(f.atom(0).position, inside);
        assert_eq!(f.atom(1).position, at_zero);
        // 越界原子：x、y 折回，z 非周期保持 1.6
        let g = cell.cartesian_to_fractional(f.atom(2).position).unwrap();
        for (k, w) in [0.25, 0.6, 1.6].iter().enumerate() {
            assert!((g[k] - w).abs() < 1e-12, "第 {k} 轴分数坐标 {}，应为 {w}", g[k]);
        }
        // 位移恰为整数个晶格矢量
        let d = cell.cartesian_to_fractional(f.atom(2).position - out).unwrap();
        assert!((d - Vector3::new(-1.0, 1.0, 0.0)).norm() < 1e-12, "{d:?}");
    }
}
