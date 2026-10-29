//! Linked-cell neighbour search in a 3-D periodic cell.
//!
//! Every distance-cutoff analysis here (g(r), angles, bond lifetimes, rotcorr, the
//! network classification) goes through it. It lived in `ferro-core` while the network
//! classification did; both are analysis (review A-1), and the periodic cell it
//! builds on, [`PeriodicCell`], stays in core for io and structure to measure with.

use ferro_core::error::{ChemError, Result};
use ferro_core::PeriodicCell;
use nalgebra::Vector3;

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
    cell: PeriodicCell,
    rcut: f64,
}

impl NeighborList {
    /// `Err` when `rcut` exceeds the cell's minimum-image bound: past it each pair
    /// is seen at one image only and neighbours in farther images are lost.
    /// Also `Err` for a cutoff that is not a positive finite number.
    pub fn build(positions: &[Vector3<f64>], cell: &PeriodicCell, rcut: f64) -> Result<Self> {
        if !(rcut.is_finite() && rcut > 0.0) {
            return Err(ChemError::ValidationError(format!(
                "neighbour cutoff must be a finite number > 0, got {rcut}"
            )));
        }
        let bound = cell.mic_bound();
        if rcut > bound {
            return Err(ChemError::ValidationError(format!(
                "cutoff {rcut:.3} A exceeds the minimum-image bound {bound:.3} A of the cell"
            )));
        }
        // 每轴 bin 数 = floor(d_i / rcut)：(Mᵀ)⁻¹ 第 i 行的模 = 1/d_i，对三斜胞也成立
        let mut nf = [0, 1, 2].map(|i| (1.0 / (rcut * cell.reciprocal().row(i).norm())).floor().max(1.0));
        // More bins than atoms only adds empty ones; a tiny cutoff asked for (d / rcut)^3 of
        // them. Shrinking every axis by the same factor keeps each bin at least rcut wide
        let cap = positions.len().max(27) as f64;
        let total = nf[0] * nf[1] * nf[2];
        if total > cap {
            let s = (total / cap).cbrt();
            nf = nf.map(|v| (v / s).floor().max(1.0));
        }
        let n = nf.map(|v| v as usize);
        let mut bins = vec![Vec::new(); n[0] * n[1] * n[2]];
        let frac = positions.iter().enumerate().map(|(i, &p)| {
            let f = cell.to_fractional(p);
            let f = [f.x.rem_euclid(1.0), f.y.rem_euclid(1.0), f.z.rem_euclid(1.0)];
            let b = [0, 1, 2].map(|k| ((f[k] * n[k] as f64) as usize).min(n[k] - 1));
            bins[b[0] + b[1] * n[0] + b[2] * n[0] * n[1]].push(i);
            f
        }).collect();
        Ok(Self { n, bins, frac, cell: cell.clone(), rcut })
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
                        let cart = self.cell.to_cartesian(df.map(|v| v - v.round()));
                        let d2 = cart.norm_squared();
                        // `j == i` above is the self-exclusion; a coincident atom (d = 0)
                        // stays, so the defect shows up downstream (review R4-6)
                        if d2 <= rcut2 {
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
    use ferro_core::Cell;
    use nalgebra::Matrix3;

    fn cubic(a: f64) -> Cell {
        Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap()
    }

    #[test]
    fn a_cutoff_past_the_minimum_image_bound_is_refused() {
        let pc = PeriodicCell::new(&cubic(6.0), [true; 3]).unwrap();
        assert!(NeighborList::build(&[Vector3::zeros()], &pc, 3.0).is_ok());
        assert!(NeighborList::build(&[Vector3::zeros()], &pc, 3.01).is_err());
    }

    #[test]
    fn a_cutoff_that_is_not_positive_is_refused() {
        let pc = PeriodicCell::new(&cubic(6.0), [true; 3]).unwrap();
        for rcut in [0.0, -1.0, f64::NAN] {
            assert!(NeighborList::build(&[Vector3::zeros()], &pc, rcut).is_err(), "rcut {rcut}");
        }
    }

    /// A tiny cutoff used to ask for (d / rcut)^3 bins: 8e15 here. The count is capped
    /// by the atom count (at least 27), and the wider bins still find the pair.
    #[test]
    fn a_tiny_cutoff_does_not_blow_up_the_bin_count() {
        let pc = PeriodicCell::new(&cubic(20.0), [true; 3]).unwrap();
        let pos = [Vector3::new(1.0, 1.0, 1.0), Vector3::new(1.0 + 5e-5, 1.0, 1.0), Vector3::new(9.0, 9.0, 9.0)];
        let nl = NeighborList::build(&pos, &pc, 1e-4).unwrap();
        assert!(nl.bins.len() <= 27, "{} bins", nl.bins.len());
        assert_eq!(nl.neighbors_of(0).iter().map(|n| n.0).collect::<Vec<_>>(), vec![1]);
        assert!(nl.neighbors_of(2).is_empty());
    }

    /// Review R4-6: a coincident atom is a neighbour at distance 0; only `i` itself is left out.
    #[test]
    fn a_coincident_atom_is_kept() {
        let pc = PeriodicCell::new(&cubic(10.0), [true; 3]).unwrap();
        let pos = [Vector3::new(1.0, 1.0, 1.0), Vector3::new(1.0, 1.0, 1.0)];
        let nl = NeighborList::build(&pos, &pc, 2.0).unwrap();
        let n = nl.neighbors_of(0);
        assert_eq!(n.len(), 1);
        assert_eq!((n[0].0, n[0].1), (1, 0.0));
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
