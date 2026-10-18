//! Periodic geometry of a 3-D periodic frame: minimum image and neighbour search.
//!
//! Every analysis that measures distances in a box goes through [`PeriodicCell`].
//! Building one is where the box is validated — singular, non-finite and not
//! fully periodic cells are refused there — so the hot loops that use it cannot
//! fail. Only 3-D periodic systems are supported: a frame that carries a box but
//! whose `pbc` is not `[true; 3]` (an LAMMPS `f`/`s` boundary, CP2K `PERIODIC XY`,
//! a DeePMD `nopbc` system) is an error, never quietly treated as a bulk.

use nalgebra::{Matrix3, Vector3};

use crate::cell::Cell;
use crate::error::{ChemError, Result};
use crate::frame::Frame;

/// A validated 3-D periodic cell with its inverse computed once.
#[derive(Debug, Clone)]
pub struct PeriodicCell {
    /// `Mᵀ`: fractional → Cartesian (`cart = Mᵀ·frac`, rows of `M` are lattice vectors)
    mt: Matrix3<f64>,
    /// `(Mᵀ)⁻¹`: Cartesian → fractional; its rows are the reciprocal lattice vectors
    inv: Matrix3<f64>,
}

impl PeriodicCell {
    /// Validates `cell` and `pbc`: the matrix must be finite and invertible and
    /// all three directions periodic.
    pub fn new(cell: &Cell, pbc: [bool; 3]) -> Result<Self> {
        if pbc != [true; 3] {
            let flag = |p: bool| if p { 'T' } else { 'F' };
            return Err(ChemError::ValidationError(format!(
                "only 3-D periodic cells are supported, got pbc = \"{} {} {}\"",
                flag(pbc[0]), flag(pbc[1]), flag(pbc[2])
            )));
        }
        if !cell.matrix.iter().all(|v| v.is_finite()) {
            return Err(ChemError::ValidationError("cell matrix is not finite".into()));
        }
        let mt = cell.matrix.transpose();
        let inv = mt
            .try_inverse()
            .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))?;
        Ok(Self { mt, inv })
    }

    /// The frame's cell; `Err` when it has none or the cell is refused by [`Self::new`].
    pub fn of(frame: &Frame) -> Result<Self> {
        let cell = frame.cell.as_ref().ok_or_else(|| {
            ChemError::ValidationError("frame has no periodic cell".into())
        })?;
        Self::new(cell, frame.pbc)
    }

    pub fn to_fractional(&self, cart: Vector3<f64>) -> Vector3<f64> {
        self.inv * cart
    }

    pub fn to_cartesian(&self, frac: Vector3<f64>) -> Vector3<f64> {
        self.mt * frac
    }

    /// Minimum-image form of a displacement: each fractional component in [-0.5, 0.5].
    /// Exact for distances below [`Self::mic_bound`].
    pub fn minimum_image(&self, diff: Vector3<f64>) -> Vector3<f64> {
        let f = self.inv * diff;
        self.mt * f.map(|c| c - c.round())
    }

    /// `(Mᵀ)⁻¹`; row `i` is reciprocal lattice vector `i` (without 2π), of norm `1/dᵢ`.
    pub fn reciprocal(&self) -> &Matrix3<f64> {
        &self.inv
    }

    /// Interplanar spacings d₁, d₂, d₃ \[Å\]; see [`Cell::interplanar_spacings`].
    pub fn spacings(&self) -> [f64; 3] {
        [0, 1, 2].map(|i| 1.0 / self.inv.row(i).norm())
    }

    /// Largest distance the minimum image measures correctly: half the smallest
    /// interplanar spacing.
    pub fn mic_bound(&self) -> f64 {
        let d = self.spacings();
        d[0].min(d[1]).min(d[2]) * 0.5
    }
}

/// Linked-cell neighbour search within a cutoff, O(N) per frame.
///
/// `ferro-structure`'s `box_builder` keeps a cell list of its own on purpose: a cubic box
/// rebuilt every relaxation step, with force sums tied to its visiting order.
pub struct NeighborList {
    n: [usize; 3],
    /// `bins[ix + iy*nx + iz*nx*ny]` = atom indices in that bin
    bins: Vec<Vec<usize>>,
    /// Fractional coordinates folded into [0, 1)
    frac: Vec<[f64; 3]>,
    mt: Matrix3<f64>,
    rcut: f64,
}

impl NeighborList {
    /// `Err` when `rcut` exceeds the cell's minimum-image bound: past it each pair
    /// is seen at one image only and neighbours in farther images are lost.
    pub fn build(positions: &[Vector3<f64>], cell: &PeriodicCell, rcut: f64) -> Result<Self> {
        let bound = cell.mic_bound();
        if rcut > bound {
            return Err(ChemError::ValidationError(format!(
                "cutoff {rcut:.3} A exceeds the minimum-image bound {bound:.3} A of the cell"
            )));
        }
        // 每轴 bin 数 = floor(d_i / rcut)：(Mᵀ)⁻¹ 第 i 行的模 = 1/d_i，对三斜胞也成立
        let n = [0, 1, 2].map(|i| ((1.0 / (rcut * cell.inv.row(i).norm())).floor() as usize).max(1));
        let mut bins = vec![Vec::new(); n[0] * n[1] * n[2]];
        let frac = positions.iter().enumerate().map(|(i, &p)| {
            let f = cell.inv * p;
            let f = [f.x.rem_euclid(1.0), f.y.rem_euclid(1.0), f.z.rem_euclid(1.0)];
            let b = [0, 1, 2].map(|k| ((f[k] * n[k] as f64) as usize).min(n[k] - 1));
            bins[b[0] + b[1] * n[0] + b[2] * n[0] * n[1]].push(i);
            f
        }).collect();
        Ok(Self { n, bins, frac, mt: cell.mt, rcut })
    }

    /// Atoms other than `i` within the cutoff, sorted by index:
    /// `(index, distance, minimum-image vector i → index)`.
    ///
    /// A superset: the test is inclusive and allows a relative `1e-9` past the cutoff,
    /// so a caller that re-measures a pair with its own arithmetic (and its own `<` or
    /// `<=`) never loses one that sits on the cutoff to a last-bit difference.
    pub fn neighbors_of(&self, i: usize) -> Vec<(usize, f64, [f64; 3])> {
        let fi = self.frac[i];
        let rcut2 = self.rcut * self.rcut * (1.0 + 2e-9);
        let [nx, ny, nz] = self.n;
        let c = [0, 1, 2].map(|k| (fi[k] * self.n[k] as f64) as i64);

        let mut out = Vec::new();
        // 每轴只有 2 个 bin 时 -1 与 +1 落在同一个 bin，去重防止重复计入
        let mut seen: Vec<usize> = Vec::with_capacity(27);
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let ix = (c[0] + dx).rem_euclid(nx as i64) as usize;
                    let iy = (c[1] + dy).rem_euclid(ny as i64) as usize;
                    let iz = (c[2] + dz).rem_euclid(nz as i64) as usize;
                    let flat = ix + iy * nx + iz * nx * ny;
                    if seen.contains(&flat) { continue; }
                    seen.push(flat);
                    for &j in &self.bins[flat] {
                        if j == i { continue; }
                        let fj = self.frac[j];
                        let df = Vector3::new(fj[0] - fi[0], fj[1] - fi[1], fj[2] - fi[2]);
                        let cart = self.mt * df.map(|v| v - v.round());
                        let d2 = cart.norm_squared();
                        if d2 > 0.0 && d2 <= rcut2 {
                            out.push((j, d2.sqrt(), [cart.x, cart.y, cart.z]));
                        }
                    }
                }
            }
        }
        out.sort_unstable_by_key(|&(j, _, _)| j);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic(a: f64) -> Cell {
        Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap()
    }

    #[test]
    fn only_finite_invertible_3d_periodic_cells_are_accepted() {
        assert!(PeriodicCell::new(&cubic(10.0), [true; 3]).is_ok());
        let slab = PeriodicCell::new(&cubic(10.0), [true, true, false]).unwrap_err();
        assert!(slab.to_string().contains("pbc = \"T T F\""), "{slab}");
        assert!(PeriodicCell::new(&cubic(10.0), [false; 3]).is_err());

        let mut flat = cubic(10.0);
        flat.matrix[(2, 2)] = 0.0;
        assert!(PeriodicCell::new(&flat, [true; 3]).unwrap_err().to_string().contains("singular"));
        let mut nan = cubic(10.0);
        nan.matrix[(0, 1)] = f64::NAN;
        assert!(PeriodicCell::new(&nan, [true; 3]).unwrap_err().to_string().contains("not finite"));

        let mut frame = Frame::new();
        assert!(PeriodicCell::of(&frame).is_err());
        frame.cell = Some(cubic(10.0));
        frame.pbc = [true; 3];
        assert!(PeriodicCell::of(&frame).is_ok());
    }

    #[test]
    fn minimum_image_matches_cell() {
        let cell = Cell::from_lengths_angles(7.0, 9.0, 11.0, 75.0, 85.0, 95.0).unwrap();
        let pc = PeriodicCell::new(&cell, [true; 3]).unwrap();
        let d = Vector3::new(6.1, -8.2, 4.4);
        assert_eq!(pc.minimum_image(d), cell.minimum_image(d).unwrap());
        assert_eq!(pc.mic_bound(), cell.minimum_image_cutoff().unwrap());
        let f = Vector3::new(0.3, 0.5, 0.7);
        assert!((pc.to_fractional(pc.to_cartesian(f)) - f).norm() < 1e-12);
    }

    #[test]
    fn a_cutoff_past_the_minimum_image_bound_is_refused() {
        let pc = PeriodicCell::new(&cubic(6.0), [true; 3]).unwrap();
        assert!(NeighborList::build(&[Vector3::zeros()], &pc, 3.0).is_ok());
        assert!(NeighborList::build(&[Vector3::zeros()], &pc, 3.01).is_err());
    }

    /// Reference pairs from ASE 3.26 `neighbor_list('ijd', atoms, 3.5)` on the same
    /// strongly triclinic cell (bound 3.78 A), 14 atoms drawn uniformly in fractional
    /// [-0.3, 1.3) so some sit outside the box.
    #[test]
    fn neighbors_match_ase_on_a_triclinic_cell() {
        let m = Matrix3::new(
            9.0, 0.0, 0.0,
            4.0, 8.5, 0.0,
            -3.0, 2.5, 9.5,
        );
        let pc = PeriodicCell::new(&Cell::from_matrix(m), [true; 3]).unwrap();
        let pos: Vec<Vector3<f64>> = TRICLINIC_POS.iter().map(|p| Vector3::new(p[0], p[1], p[2])).collect();
        let nl = NeighborList::build(&pos, &pc, 3.5).unwrap();
        let mut got: Vec<(usize, usize, f64)> = (0..pos.len())
            .flat_map(|i| nl.neighbors_of(i).into_iter().map(move |(j, d, _)| (i, j, d)))
            .collect();
        got.sort_by_key(|g| (g.0, g.1));
        assert_eq!(got.len(), TRICLINIC_PAIRS.len());
        for (g, e) in got.iter().zip(TRICLINIC_PAIRS) {
            assert_eq!((g.0, g.1), (e.0, e.1));
            assert!((g.2 - e.2).abs() < 1e-9, "{g:?} vs {e:?}");
        }
    }

    const TRICLINIC_POS: [[f64; 3]; 14] = [
        [8.02025173213556, 12.004850454167002, 8.940422491726943],
        [-2.029008778605691, 4.276475256377713, 10.42801237002318],
        [-1.494250994615842, 11.056984205013807, 9.265455317031105],
        [4.341227914509473, 1.9349434531457577, 1.3820693039317549],
        [1.096778776870157, 4.771230795835807, 4.819133536160891],
        [7.536786471626791, 13.409451531562754, 9.198461172049047],
        [11.255244097585932, 11.011092801416032, 0.42269221318110434],
        [3.016385116685719, 5.206306649958753, -2.182081478986973],
        [-1.4287044574848329, 5.567312056991793, 4.236331584944395],
        [11.766499259041876, 7.313947647475798, 4.964588228312611],
        [5.682461648037352, 0.11337904174172311, -2.670730811753911],
        [3.2366841718946673, 6.914063739940994, 0.19922220460232723],
        [-1.6390070812433164, 0.0709766111152157, 9.766725492986534],
        [-3.2887192226150166, 3.860679157990797, 10.531048740508595],
    ];
    const TRICLINIC_PAIRS: [(usize, usize, f64); 60] = [
        (0, 1, 3.3934041926754204),
        (0, 2, 1.1264145860399621),
        (0, 3, 3.436893065764812),
        (0, 5, 1.5077223899917684),
        (0, 6, 1.8135162602391695),
        (0, 10, 2.653558339113941),
        (0, 11, 2.961023587100929),
        (0, 13, 2.3486773217398955),
        (1, 0, 3.3934041926754204),
        (1, 3, 3.4043758952614516),
        (1, 6, 2.8583725857316957),
        (1, 13, 1.3305537538268648),
        (2, 0, 1.1264145860399621),
        (2, 3, 3.2901859025180435),
        (2, 5, 2.3536257268845846),
        (2, 6, 2.6488397901403413),
        (2, 10, 2.5721557698672095),
        (2, 11, 2.4255881591834703),
        (2, 13, 2.857571983567912),
        (3, 0, 3.436893065764812),
        (3, 1, 3.4043758952614516),
        (3, 2, 3.2901859025180435),
        (3, 5, 2.807164340745951),
        (3, 6, 3.1215129006622626),
        (4, 8, 2.7113592634901043),
        (4, 9, 3.0454117198938304),
        (4, 10, 3.3483368558650035),
        (5, 0, 1.5077223899917684),
        (5, 2, 2.3536257268845846),
        (5, 3, 2.807164340745951),
        (5, 6, 1.0251940835644024),
        (5, 10, 3.4080625638081044),
        (5, 13, 2.757560263159354),
        (6, 0, 1.8135162602391695),
        (6, 1, 2.8583725857316957),
        (6, 2, 2.6488397901403413),
        (6, 3, 3.1215129006622626),
        (6, 5, 1.0251940835644024),
        (6, 13, 1.9528420387325582),
        (7, 11, 2.9386345784371586),
        (7, 12, 3.498784856667299),
        (8, 4, 2.7113592634901043),
        (9, 4, 3.0454117198938304),
        (10, 0, 2.653558339113941),
        (10, 2, 2.5721557698672095),
        (10, 4, 3.3483368558650035),
        (10, 5, 3.4080625638081044),
        (11, 0, 2.961023587100929),
        (11, 2, 2.4255881591834703),
        (11, 7, 2.9386345784371586),
        (11, 12, 2.286490816235405),
        (11, 13, 3.0983414507877547),
        (12, 7, 3.498784856667299),
        (12, 11, 2.286490816235405),
        (13, 0, 2.3486773217398955),
        (13, 1, 1.3305537538268648),
        (13, 2, 2.857571983567912),
        (13, 5, 2.757560263159354),
        (13, 6, 1.9528420387325582),
        (13, 11, 3.0983414507877547),
    ];
}
