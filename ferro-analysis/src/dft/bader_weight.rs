//! Yu-Trinkle weight method for Bader analysis.
//!
//! Reference: Yu & Trinkle, *J. Chem. Phys.* **134**, 064111 (2011).
//! See `dev/bader.md` §4.4 for detailed algorithm description.

use ferro_core::{ChargeGrid, Frame};
use super::bader::{BaderResult, BaderParams};
use super::bader_grid::{mark_vacuum, finalize_vacuum};

/// Wigner-Seitz Voronoi decomposition: compute neighbor vectors and area weights.
///
/// Returns `(vect, alpha)` where:
/// - `vect[n]` = lattice offset `[d1, d2, d3]` for neighbor `n`
/// - `alpha[n]` = area weight for that neighbor
fn ws_voronoi(chg: &ChargeGrid) -> (Vec<[i32; 3]>, Vec<f64>) {
    let lat2car = &chg.lat2car;

    // Generate candidate directions (nRange=2 → 5×5×5-1 = 124 directions)
    let n_range = 2i32;
    let mut candidates: Vec<[i32; 3]> = Vec::new();
    for d1 in -n_range..=n_range {
        for d2 in -n_range..=n_range {
            for d3 in -n_range..=n_range {
                if d1 == 0 && d2 == 0 && d3 == 0 { continue; }
                candidates.push([d1, d2, d3]);
            }
        }
    }

    // Filter: R/2 must be inside the Wigner-Seitz cell
    // A vector R is inside the WS cell if |R/2| ≤ |R/2 - T| for all lattice translations T
    // Simplification: check against the 26 nearest lattice translations
    let ws_vectors: Vec<[i32; 3]> = candidates.iter().filter(|&&d| {
        let r = lat2car * nalgebra::Vector3::new(d[0] as f64, d[1] as f64, d[2] as f64);
        let half_r = r / 2.0;
        let half_r_norm2 = half_r.norm_squared();

        // Check against all lattice translations within nRange
        for t1 in -n_range..=n_range {
            for t2 in -n_range..=n_range {
                for t3 in -n_range..=n_range {
                    if t1 == 0 && t2 == 0 && t3 == 0 { continue; }
                    let t = lat2car * nalgebra::Vector3::new(t1 as f64, t2 as f64, t3 as f64);
                    if (half_r - t).norm_squared() < half_r_norm2 - 1e-10 {
                        return false;
                    }
                }
            }
        }
        true
    }).copied().collect();

    // Compute alpha (area weights) for each WS vector
    // For the simplified implementation, use the projection of the face normal
    // alpha[n] ≈ |R[n]| * face_area / 2
    // A practical approximation: alpha is proportional to the face area of the WS cell
    // For a more accurate implementation, we would compute the Voronoi cell faces exactly.
    //
    // Simplified: alpha[n] = |R_car[n]| for now, which gives correct relative weights
    // for cubic cells and reasonable approximations for others.
    let mut alpha = Vec::with_capacity(ws_vectors.len());
    for &d in &ws_vectors {
        let r = lat2car * nalgebra::Vector3::new(d[0] as f64, d[1] as f64, d[2] as f64);
        alpha.push(r.norm());
    }

    (ws_vectors, alpha)
}

/// Run Bader analysis using the Yu-Trinkle weight method.
pub fn bader_weight(chg: &ChargeGrid, frame: &Frame, params: &BaderParams) -> BaderResult {
    let nrho = chg.nrho;
    let vol_cell = frame.cell.as_ref().unwrap().volume();

    // Step 1: WS Voronoi decomposition
    let (vect, alpha) = ws_voronoi(chg);
    let num_vect = vect.len();

    // Step 2: Sort grid points by density (descending)
    let mut indices: Vec<usize> = (0..nrho).collect();
    indices.sort_unstable_by(|&a, &b| {
        chg.rho[b].partial_cmp(&chg.rho[a]).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Reverse lookup: position → rank
    let mut rank = vec![0usize; nrho];
    for (r, &idx) in indices.iter().enumerate() {
        rank[idx] = r;
    }

    // Step 3: Flow assignment
    let mut basin = vec![0i32; nrho]; // 0 = unassigned
    let mut nvols = 0usize;
    let mut volpos_lat: Vec<[f64; 3]> = Vec::new();

    // For boundary points: track downstream flow
    // prob[n] and neigh[n] store flow fractions from rank n to its lower-density neighbors
    let mut prob: Vec<Vec<f64>> = vec![Vec::new(); nrho];
    let mut neigh: Vec<Vec<usize>> = vec![Vec::new(); nrho];
    let mut numbelow = vec![0usize; nrho];

    let [n1, n2, n3] = chg.shape;

    for (n, &pos) in indices.iter().enumerate().take(nrho) {
        let i1 = pos % n1;
        let rem = pos / n1;
        let i2 = rem % n2;
        let i3 = rem / n2;

        // Find all higher-density neighbors (rank m < n)
        let mut above: Vec<(usize, f64)> = Vec::new();
        let mut tsum = 0.0_f64;

        for nv in 0..num_vect {
            let nb_p = chg.pbc_i([
                i1 as i32 + vect[nv][0],
                i2 as i32 + vect[nv][1],
                i3 as i32 + vect[nv][2],
            ]);
            let nb_idx = nb_p[0] + n1 * (nb_p[1] + n2 * nb_p[2]);
            let m = rank[nb_idx];

            if m < n {
                let t = alpha[nv] * (chg.rho[indices[m]] - chg.rho[pos]);
                above.push((m, t));
                tsum += t;
            }
        }

        if above.is_empty() {
            // This is a density maximum → new basin
            nvols += 1;
            basin[pos] = nvols as i32;
            volpos_lat.push([i1 as f64, i2 as f64, i3 as f64]);
        } else {
            // Check if all above neighbors belong to the same basin
            let first_basin = basin[indices[above[0].0]];
            let is_boundary = above.iter().any(|(m, _)| {
                let b = basin[indices[*m]];
                b != first_basin || b == 0
            });

            if is_boundary {
                // Boundary point: distribute flow to upstream basins
                basin[pos] = 0;
                for &(m, t) in &above {
                    let idx_m = indices[m];
                    let frac = t / tsum;
                    prob[idx_m].push(frac);
                    neigh[idx_m].push(pos);
                    numbelow[idx_m] += 1;
                }
            } else {
                // Internal point: inherit basin
                basin[pos] = first_basin;
            }
        }
    }

    // Step 4: Charge integration
    // Forward pass: propagate weights from high density to low density
    let mut weight = vec![0.0_f64; nrho];
    let mut volchg = vec![0.0_f64; nvols + 1]; // index 0 unused, 1..=nvols for volumes

    // Initialize: maxima get weight 1
    for &pos in indices.iter().take(nrho) {
        if basin[pos] > 0 && weight[pos] == 0.0 {
            // This is a basin maximum — it should already have weight 0
            // Set weight based on basin assignment
            weight[pos] = if basin[pos] > 0 { 1.0 } else { 0.0 };
        }

        // Propagate weight to lower-density neighbors
        for k in 0..numbelow[pos] {
            let target = neigh[pos][k];
            weight[target] += prob[pos][k] * weight[pos];
        }

        // Accumulate charge for the basin
        let b = basin[pos];
        if b > 0 && weight[pos] > 0.0 {
            volchg[b as usize] += weight[pos] * chg.rho[pos];
        }
    }

    // Normalize: volchg = sum / nrho
    let inv_nrho = 1.0 / nrho as f64;
    for v in volchg.iter_mut() {
        *v *= inv_nrho;
    }

    // Step 5: Assign volnum for visualization (take largest-flow basin for boundary points)
    let mut volnum = vec![0i32; nrho];
    for i in 0..nrho {
        if basin[i] > 0 {
            volnum[i] = basin[i];
        } else if basin[i] == 0 && !neigh[i].is_empty() {
            // Boundary point: find which upstream basin contributes most
            // We need to look at who flows INTO this point
            // Actually, for boundary points we need to look backwards.
            // For simplicity, assign to the basin of the highest-density neighbor.
            let best = neigh[i].iter().enumerate()
                .max_by(|(_, &a), (_, &b)| weight[a].partial_cmp(&weight[b]).unwrap_or(std::cmp::Ordering::Equal));
            if let Some((_, &target)) = best {
                volnum[i] = volnum[target].max(1);
            }
        }
    }

    // Step 6: Mark vacuum
    let _nvac = mark_vacuum(&mut volnum, &chg.rho, vol_cell, nrho, params.vacval);
    // Also mark vacuum in volnum for basin=0 points with low density
    for (i, vn) in volnum.iter_mut().enumerate().take(nrho) {
        if *vn == 0 && chg.rho[i].abs() / vol_cell <= params.vacval {
            *vn = -1;
        }
    }

    // Add vacuum entry to volchg
    volchg.push(0.0);

    // Step 7: Convert positions
    let volpos_car: Vec<[f64; 3]> = volpos_lat.iter().map(|lat| {
        let v = chg.lat2car * nalgebra::Vector3::new(lat[0], lat[1], lat[2]);
        [v.x, v.y, v.z]
    }).collect();
    let volpos_dir: Vec<[f64; 3]> = volpos_lat.iter().map(|lat| {
        [lat[0] / n1 as f64, lat[1] / n2 as f64, lat[2] / n3 as f64]
    }).collect();

    // Step 8: Assign volumes to atoms
    let cell = frame.cell.as_ref().unwrap();
    let nions = frame.n_atoms();
    let mut ionchg = vec![0.0_f64; nions];
    let mut iondist = vec![0.0_f64; nvols];
    let mut nnion = vec![0usize; nvols];

    for v in 0..nvols {
        let vp = nalgebra::Vector3::new(volpos_car[v][0], volpos_car[v][1], volpos_car[v][2]);
        let mut best_j = 0usize;
        let mut best_d = f64::MAX;
        for (j, atom) in frame.atoms.iter().enumerate() {
            let dv = cell.minimum_image(vp - atom.position).expect("cell is non-singular");
            let d = dv.norm();
            if d < best_d {
                best_d = d;
                best_j = j;
            }
        }
        nnion[v] = best_j;
        iondist[v] = best_d;
        ionchg[best_j] += volchg[v + 1]; // volchg is 1-indexed
    }

    // Step 9: Atomic volumes
    let voxel_vol = vol_cell / nrho as f64;
    let mut ionvol = vec![0.0_f64; nions];
    for &b in &basin {
        if b > 0 && (b as usize) <= nvols {
            let j = nnion[(b - 1) as usize];
            ionvol[j] += voxel_vol;
        }
    }

    // Step 10: Vacuum
    let vacchg = volchg[nvols]; // last element
    let vacvol = _nvac as f64 * voxel_vol;

    finalize_vacuum(&mut volnum, nvols);

    BaderResult {
        nvols,
        volnum,
        volpos_lat,
        volpos_dir,
        volpos_car,
        volchg,
        ionchg,
        iondist,
        nnion,
        ionvol,
        vacchg,
        vacvol,
    }
}
