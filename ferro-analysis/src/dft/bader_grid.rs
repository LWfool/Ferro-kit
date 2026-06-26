//! Grid-based Bader methods: on-grid, near-grid, off-grid gradient ascent.
//!
//! Reference: Henkelman et al., *Comput. Mater. Sci.* **36**, 254 (2006).

use ferro_core::{Cell, ChargeGrid, Frame};
use nalgebra::Vector3;
use super::bader::{BaderResult, BaderParams};

// ─── Shared utilities ────────────────────────────────────────────────────────

/// Assign each Bader volume to its nearest atom (considering PBC).
fn assign_chg2atom(
    volpos_car: &[[f64; 3]],
    volchg: &[f64],
    frame: &Frame,
    cell: &Cell,
    nvols: usize,
) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let nions = frame.n_atoms();
    let mut ionchg = vec![0.0_f64; nions];
    let mut iondist = vec![0.0_f64; nvols];
    let mut nnion = vec![0_usize; nvols];

    for v in 0..nvols {
        let vp = Vector3::new(volpos_car[v][0], volpos_car[v][1], volpos_car[v][2]);
        let mut best_j = 0usize;
        let mut best_d = f64::MAX;
        for (j, atom) in frame.atoms.iter().enumerate() {
            let dv = cell.minimum_image(vp - atom.position)
                .expect("cell is non-singular");
            let d = dv.norm();
            if d < best_d {
                best_d = d;
                best_j = j;
            }
        }
        nnion[v] = best_j;
        iondist[v] = best_d;
        ionchg[best_j] += volchg[v];
    }
    (ionchg, iondist, nnion)
}

/// Compute Bader volume for each atom (Å³).
fn calc_atomic_vol(volnum: &[i32], nvols: usize, nnion: &[usize], nions: usize, vol_cell: f64, nrho: usize) -> Vec<f64> {
    let mut ionvol = vec![0.0_f64; nions];
    let voxel_vol = vol_cell / nrho as f64;
    for &vn in volnum.iter() {
        if vn > 0 && (vn as usize) <= nvols {
            let j = nnion[(vn - 1) as usize];
            ionvol[j] += voxel_vol;
        }
    }
    ionvol
}

/// Compute volume charges: `volchg[v] = Σ rho_stored / nrho` for each Bader volume.
fn calc_volchg(volnum: &[i32], rho: &[f64], nvols: usize, nrho: usize) -> Vec<f64> {
    let mut volchg = vec![0.0_f64; nvols];
    for (i, &vn) in volnum.iter().enumerate() {
        if vn > 0 && (vn as usize) <= nvols {
            volchg[(vn - 1) as usize] += rho[i];
        }
    }
    let inv_nrho = 1.0 / nrho as f64;
    for v in volchg.iter_mut() {
        *v *= inv_nrho;
    }
    volchg
}

/// Convert volpos lattice coordinates to Cartesian and fractional.
fn convert_volpos(
    volpos_lat: &[[f64; 3]],
    chg: &ChargeGrid,
    nvols: usize,
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut volpos_car = vec![[0.0_f64; 3]; nvols];
    let mut volpos_dir = vec![[0.0_f64; 3]; nvols];
    for v in 0..nvols {
        let lat = Vector3::new(volpos_lat[v][0], volpos_lat[v][1], volpos_lat[v][2]);
        let car = chg.lat2car * lat;
        volpos_car[v] = [car.x, car.y, car.z];
        // fractional = lat / N_i
        let dir = Vector3::new(
            lat.x / chg.shape[0] as f64,
            lat.y / chg.shape[1] as f64,
            lat.z / chg.shape[2] as f64,
        );
        volpos_dir[v] = [dir.x, dir.y, dir.z];
    }
    (volpos_car, volpos_dir)
}

/// Mark vacuum grid points: `|rho_stored| / V_cell <= vacval`.
pub(crate) fn mark_vacuum(volnum: &mut [i32], rho: &[f64], vol_cell: f64, _nrho: usize, vacval: f64) -> usize {
    let mut nvac = 0usize;
    for (i, vn) in volnum.iter_mut().enumerate() {
        if rho[i].abs() / vol_cell <= vacval {
            *vn = -1;
            nvac += 1;
        }
    }
    nvac
}

/// Finalize vacuum: replace -1 with nvols+1.
pub(crate) fn finalize_vacuum(volnum: &mut [i32], nvols: usize) {
    let vac_id = (nvols + 1) as i32;
    for vn in volnum.iter_mut() {
        if *vn == -1 {
            *vn = vac_id;
        }
    }
}

// ─── On-grid method ──────────────────────────────────────────────────────────

/// One step of on-grid gradient ascent: move to the neighbor with maximum
/// distance-corrected density.  Returns `true` if the point moved.
fn step_ongrid(chg: &ChargeGrid, p: &mut [usize; 3]) -> bool {
    let rho_ctr = chg.rho_val([p[0] as i32, p[1] as i32, p[2] as i32]);
    let mut rho_max = rho_ctr;
    let mut pm = *p;

    for d1 in -1i32..=1 {
        for d2 in -1i32..=1 {
            for d3 in -1i32..=1 {
                let pt = [p[0] as i32 + d1, p[1] as i32 + d2, p[2] as i32 + d3];
                let rho_nbr = chg.rho_val(pt);
                let li = chg.lat_i_dist_i(d1, d2, d3);
                // Distance-corrected extrapolated density at unit distance
                let rho_tmp = rho_ctr + (rho_nbr - rho_ctr) * li;
                if rho_tmp > rho_max {
                    rho_max = rho_tmp;
                    pm = chg.pbc_i(pt);
                }
            }
        }
    }

    let moved = pm != *p;
    *p = pm;
    moved
}

/// Climb from `start` to a local density maximum using on-grid steps.
fn max_ongrid(chg: &ChargeGrid, start: [usize; 3]) -> [usize; 3] {
    let mut p = start;
    // Safety limit to prevent infinite loops on pathological grids
    let max_iter = chg.nrho;
    for _ in 0..max_iter {
        if !step_ongrid(chg, &mut p) {
            break;
        }
    }
    p
}

/// Run Bader analysis using the on-grid method.
pub fn bader_ongrid(chg: &ChargeGrid, frame: &Frame, params: &BaderParams) -> BaderResult {
    let [n1, n2, n3] = chg.shape;
    let nrho = chg.nrho;
    let vol_cell = frame.cell.as_ref().unwrap().volume();

    let mut volnum = vec![0i32; nrho];
    let mut nvols = 0usize;
    let mut volpos_lat: Vec<[f64; 3]> = Vec::new();

    // 1. Gradient ascent from every unassigned grid point, tracing the full
    //    path so every point along it is assigned to the maximum's volume.
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let idx = i1 + n1 * (i2 + n2 * i3);
                if volnum[idx] != 0 {
                    continue;
                }
                // Trace the full path
                let mut path: Vec<usize> = Vec::new();
                let mut p = [i1, i2, i3];
                let max_iter = chg.nrho;
                for _ in 0..max_iter {
                    let pi = p[0] + n1 * (p[1] + n2 * p[2]);
                    path.push(pi);
                    if !step_ongrid(chg, &mut p) {
                        break;
                    }
                }
                let idx_max = p[0] + n1 * (p[1] + n2 * p[2]);

                let vol_id = if volnum[idx_max] > 0 {
                    volnum[idx_max]
                } else {
                    nvols += 1;
                    volpos_lat.push([p[0] as f64, p[1] as f64, p[2] as f64]);
                    nvols as i32
                };

                for &pi in &path {
                    if volnum[pi] == 0 {
                        volnum[pi] = vol_id;
                    }
                }
                if volnum[idx_max] == 0 {
                    volnum[idx_max] = vol_id;
                }
            }
        }
    }

    // 2. Mark vacuum points
    let _nvac = mark_vacuum(&mut volnum, &chg.rho, vol_cell, nrho, params.vacval);

    // 3. Compute volume charges
    let mut volchg = calc_volchg(&volnum, &chg.rho, nvols, nrho);
    // Zero out vacuum volume charge (index nvols, i.e. nvols+1-th entry)
    volchg.push(0.0);

    // 4. Convert volume positions
    let (volpos_car, volpos_dir) = convert_volpos(&volpos_lat, chg, nvols);

    // 5. Assign volumes to atoms
    let cell = frame.cell.as_ref().unwrap();
    let (ionchg, iondist, nnion) = assign_chg2atom(&volpos_car, &volchg, frame, cell, nvols);

    // 6. Atomic volumes
    let ionvol = calc_atomic_vol(&volnum, nvols, &nnion, frame.n_atoms(), vol_cell, nrho);

    // 7. Vacuum charge/volume
    let vac_vol = nvols; // index in volchg for vacuum
    let vacchg = if vac_vol < volchg.len() { volchg[vac_vol] } else { 0.0 };
    let vacvol = _nvac as f64 * vol_cell / nrho as f64;

    // 8. Finalize vacuum numbering
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

// ─── Near-grid method ────────────────────────────────────────────────────────

/// Compute gradient direction using central differences with maximum plateau suppression.
///
/// For each axis i: `grad[i] = (rho(p+e_i) - rho(p-e_i)) / 2`.
/// If both neighbors are lower than center (local maximum plateau), set `grad[i] = 0`.
fn rho_grad_dir(chg: &ChargeGrid, p: [usize; 3]) -> Vector3<f64> {
    let pi = [p[0] as i32, p[1] as i32, p[2] as i32];
    let rho_ctr = chg.rho_val(pi);
    let mut grad = Vector3::zeros();

    for i in 0..3 {
        let mut ei = [0i32; 3];
        ei[i] = 1;
        let rho_plus = chg.rho_val([pi[0] + ei[0], pi[1] + ei[1], pi[2] + ei[2]]);
        let rho_minus = chg.rho_val([pi[0] - ei[0], pi[1] - ei[1], pi[2] - ei[2]]);

        // If both sides are lower → local maximum plateau, suppress this component
        if rho_plus < rho_ctr && rho_minus < rho_ctr {
            grad[i] = 0.0;
        } else {
            grad[i] = (rho_plus - rho_minus) / 2.0;
        }
    }

    // Coordinate transform: lattice gradient → Cartesian → direction in lattice space
    // Fortran: rho_grad_car = grad_lat @ car2lat; rho_grad_dir = car2lat @ rho_grad_car
    // This is an approximation for non-orthogonal cells (see bader.md §4.5)
    let grad_car = chg.car2lat.transpose() * grad;
    chg.car2lat * grad_car
}

/// One step of near-grid gradient ascent with persistent `dr` state.
///
/// Returns `true` if the point moved.  Modifies `p` in place and accumulates
/// fractional displacement in `dr` (must be zeroed at path start).
fn step_neargrid(chg: &ChargeGrid, p: &mut [usize; 3], dr: &mut Vector3<f64>, known: &mut [u8]) -> bool {
    let grad = rho_grad_dir(chg, *p);
    let max_grad = grad.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);

    if max_grad < 1e-30 {
        // Gradient is essentially zero — check if this is a true maximum
        // (all 26 neighbors have lower or equal density)
        let rho_ctr = chg.rho_val([p[0] as i32, p[1] as i32, p[2] as i32]);
        let is_max = {
            let mut dominated = true;
            for d1 in -1i32..=1 {
                for d2 in -1i32..=1 {
                    for d3 in -1i32..=1 {
                        if d1 == 0 && d2 == 0 && d3 == 0 { continue; }
                        if chg.rho_val([p[0] as i32 + d1, p[1] as i32 + d2, p[2] as i32 + d3]) > rho_ctr {
                            dominated = false;
                            break;
                        }
                    }
                    if !dominated { break; }
                }
                if !dominated { break; }
            }
            dominated
        };
        if is_max {
            *dr = Vector3::zeros();
            return false;
        }
        // Not a maximum — fallback to on-grid step
        let moved = step_ongrid(chg, p);
        if !moved {
            return false;
        }
        *dr = Vector3::zeros();
        return true;
    }

    // Normalize so the largest component = 1
    let coeff = 1.0 / max_grad;
    let gradrl = grad * coeff;

    // Integer step
    let rounded = Vector3::new(gradrl.x.round(), gradrl.y.round(), gradrl.z.round());
    let pm = [
        (p[0] as i32 + rounded.x as i32),
        (p[1] as i32 + rounded.y as i32),
        (p[2] as i32 + rounded.z as i32),
    ];

    // Accumulate fractional part
    *dr += gradrl - rounded;
    // Extra integer step if accumulated fraction >= 0.5
    let dr_round = Vector3::new(dr.x.round(), dr.y.round(), dr.z.round());
    let pm = [
        pm[0] + dr_round.x as i32,
        pm[1] + dr_round.y as i32,
        pm[2] + dr_round.z as i32,
    ];
    *dr -= dr_round;

    // Mark current point on path
    let pidx = p[0] + chg.shape[0] * (p[1] + chg.shape[1] * p[2]);
    known[pidx] = 1;

    // PBC fold
    let pm_usize = chg.pbc_i(pm);

    // Cycle detection: if target is already on current path, fallback to on-grid
    let pm_idx = pm_usize[0] + chg.shape[0] * (pm_usize[1] + chg.shape[1] * pm_usize[2]);
    if known[pm_idx] == 1 {
        let moved = step_ongrid(chg, p);
        *dr = Vector3::zeros();
        return moved;
    }

    let moved = pm_usize != *p;
    *p = pm_usize;
    moved
}

/// Climb from `start` to a local density maximum using near-grid steps.
/// Returns `(max_position, volnum_at_max)`.  If `known[max]==2`, the volume
/// is inherited from the already-determined point.
#[allow(dead_code)]
fn max_neargrid(chg: &ChargeGrid, start: [usize; 3], known: &mut [u8]) -> ([usize; 3], i32) {
    let mut p = start;
    let mut dr = Vector3::zeros();
    let max_iter = chg.nrho;

    for _ in 0..max_iter {
        // Check if current point is already fully determined
        let pidx = p[0] + chg.shape[0] * (p[1] + chg.shape[1] * p[2]);
        if known[pidx] == 2 {
            // Inherit volume from this known point
            return (p, -2); // sentinel: caller must look up volnum
        }

        let moved = step_neargrid(chg, &mut p, &mut dr, known);
        if !moved {
            return (p, 0); // new maximum
        }
    }
    (p, 0)
}

/// Mark points whose 6 face-neighbors all belong to the same volume as `known=2`.
fn known_volnum_ongrid(chg: &ChargeGrid, volnum: &[i32], known: &mut [u8]) {
    let [n1, n2, n3] = chg.shape;
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let idx = i1 + n1 * (i2 + n2 * i3);
                if known[idx] == 2 { continue; }
                let vn = volnum[idx];
                if vn <= 0 { continue; }

                let pi = [i1 as i32, i2 as i32, i3 as i32];
                let face_offsets = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
                let all_same = face_offsets.iter().all(|off| {
                    let nb = chg.pbc_i([pi[0]+off[0], pi[1]+off[1], pi[2]+off[2]]);
                    let nidx = nb[0] + n1 * (nb[1] + n2 * nb[2]);
                    volnum[nidx] == vn
                });
                if all_same {
                    known[idx] = 2;
                }
            }
        }
    }
}

/// Edge refinement: re-assign boundary points whose volume might be wrong.
fn refine_edge(chg: &ChargeGrid, volnum: &mut [i32], _known: &mut [u8], nvols: &mut usize,
               volpos_lat: &mut Vec<[f64; 3]>, refine_count: i32) {
    let [n1, n2, n3] = chg.shape;

    let mut single_pass = || -> usize {
        // 1. Mark edge points: volnum[p] = -volnum[p] if any face neighbor differs
        let mut edge_points = Vec::new();
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                for i1 in 0..n1 {
                    let idx = i1 + n1 * (i2 + n2 * i3);
                    let vn = volnum[idx];
                    if vn <= 0 { continue; }

                    let pi = [i1 as i32, i2 as i32, i3 as i32];
                    let face_offsets = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
                    let is_edge = face_offsets.iter().any(|off| {
                        let nb = chg.pbc_i([pi[0]+off[0], pi[1]+off[1], pi[2]+off[2]]);
                        let nidx = nb[0] + n1 * (nb[1] + n2 * nb[2]);
                        volnum[nidx] != vn
                    });
                    if is_edge {
                        edge_points.push(idx);
                    }
                }
            }
        }

        // 2. Negate edge points
        for &idx in &edge_points {
            volnum[idx] = -volnum[idx];
        }

        // 3. Re-assign negated points
        let mut n_reassigned = 0usize;
        for &idx in &edge_points {
            let i1 = idx % n1;
            let rem = idx / n1;
            let i2 = rem % n2;
            let i3 = rem / n2;

            let p_max = max_ongrid(chg, [i1, i2, i3]);
            let idx_max = p_max[0] + n1 * (p_max[1] + n2 * p_max[2]);

            let new_vol = if volnum[idx_max] > 0 {
                volnum[idx_max]
            } else if volnum[idx_max] < 0 && volnum[idx_max] != -1 {
                -volnum[idx_max] // un-negate
            } else {
                *nvols += 1;
                volpos_lat.push([p_max[0] as f64, p_max[1] as f64, p_max[2] as f64]);
                *nvols as i32
            };

            let old_vol = -volnum[idx]; // was negated
            if new_vol != old_vol {
                n_reassigned += 1;
            }
            volnum[idx] = new_vol;
        }

        n_reassigned
    };

    match refine_count {
        -1 => {
            // Auto: repeat until no reassignments
            loop {
                if single_pass() == 0 { break; }
            }
        }
        -2 => {
            // Single pass only
            single_pass();
        }
        n if n > 0 => {
            for _ in 0..n {
                single_pass();
            }
        }
        _ => {}
    }
}

/// Run Bader analysis using the near-grid method (default).
pub fn bader_neargrid(chg: &ChargeGrid, frame: &Frame, params: &BaderParams) -> BaderResult {
    let [n1, n2, n3] = chg.shape;
    let nrho = chg.nrho;
    let vol_cell = frame.cell.as_ref().unwrap().volume();

    let mut volnum = vec![0i32; nrho];
    let mut known = vec![0u8; nrho];
    let mut nvols = 0usize;
    let mut volpos_lat: Vec<[f64; 3]> = Vec::new();

    // 1. Gradient ascent from every unassigned grid point
    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let idx = i1 + n1 * (i2 + n2 * i3);
                if volnum[idx] != 0 {
                    continue;
                }

                // Trace the full path with near-grid steps
                let mut path: Vec<usize> = Vec::new();
                let mut p = [i1, i2, i3];
                let mut dr = Vector3::zeros();
                let max_iter = chg.nrho;

                for _ in 0..max_iter {
                    let pi = p[0] + n1 * (p[1] + n2 * p[2]);
                    path.push(pi);

                    // Check known=2 early exit
                    if known[pi] == 2 && volnum[pi] > 0 {
                        break;
                    }

                    let moved = step_neargrid(chg, &mut p, &mut dr, &mut known);
                    if !moved {
                        break;
                    }
                }

                let idx_max = p[0] + n1 * (p[1] + n2 * p[2]);
                let vol_id = if volnum[idx_max] > 0 {
                    volnum[idx_max]
                } else {
                    nvols += 1;
                    volpos_lat.push([p[0] as f64, p[1] as f64, p[2] as f64]);
                    nvols as i32
                };

                // Assign path points
                for &pi in &path {
                    if volnum[pi] == 0 {
                        volnum[pi] = vol_id;
                    }
                    known[pi] = 0; // reset path marks
                }
                if volnum[idx_max] == 0 {
                    volnum[idx_max] = vol_id;
                }
            }
        }
    }

    // 2. Mark known=2 points
    known_volnum_ongrid(chg, &volnum, &mut known);

    // 3. Edge refinement
    if params.refine != 0 {
        refine_edge(chg, &mut volnum, &mut known, &mut nvols, &mut volpos_lat, params.refine);
    }

    // 4. Mark vacuum
    let _nvac = mark_vacuum(&mut volnum, &chg.rho, vol_cell, nrho, params.vacval);

    // 5. Volume charges
    let mut volchg = calc_volchg(&volnum, &chg.rho, nvols, nrho);
    volchg.push(0.0); // vacuum

    // 6. Convert positions
    let (volpos_car, volpos_dir) = convert_volpos(&volpos_lat, chg, nvols);

    // 7. Assign to atoms
    let cell = frame.cell.as_ref().unwrap();
    let (ionchg, iondist, nnion) = assign_chg2atom(&volpos_car, &volchg, frame, cell, nvols);

    // 8. Atomic volumes
    let ionvol = calc_atomic_vol(&volnum, nvols, &nnion, frame.n_atoms(), vol_cell, nrho);

    // 9. Vacuum
    let vacchg = if nvols < volchg.len() { volchg[nvols] } else { 0.0 };
    let vacvol = _nvac as f64 * vol_cell / nrho as f64;

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

// ─── Off-grid method ─────────────────────────────────────────────────────────

/// Trilinear interpolation of density and Cartesian gradient at fractional
/// lattice position `r = [f1, f2, f3]` (continuous, not integer).
///
/// Returns `(rho_interpolated, grad_car)`.
fn rho_grad(chg: &ChargeGrid, r: [f64; 3]) -> (f64, Vector3<f64>) {
    // Integer part (floor) and fractional part
    let p0 = [r[0].floor() as i32, r[1].floor() as i32, r[2].floor() as i32];
    let f = [r[0] - p0[0] as f64, r[1] - p0[1] as f64, r[2] - p0[2] as f64];
    let g = [1.0 - f[0], 1.0 - f[1], 1.0 - f[2]];

    // 8 corner values (with PBC)
    let rho = |di: [i32; 3]| -> f64 { chg.rho_val([p0[0]+di[0], p0[1]+di[1], p0[2]+di[2]]) };
    let rho000 = rho([0,0,0]);
    let rho001 = rho([0,0,1]);
    let rho010 = rho([0,1,0]);
    let rho011 = rho([0,1,1]);
    let rho100 = rho([1,0,0]);
    let rho101 = rho([1,0,1]);
    let rho110 = rho([1,1,0]);
    let rho111 = rho([1,1,1]);

    // Interpolate along axis 3 (z)
    let rho00_ = rho000 * g[2] + rho001 * f[2];
    let rho01_ = rho010 * g[2] + rho011 * f[2];
    let rho10_ = rho100 * g[2] + rho101 * f[2];
    let rho11_ = rho110 * g[2] + rho111 * f[2];

    // Interpolate along axis 2 (y)
    let rho0__ = rho00_ * g[1] + rho01_ * f[1];
    let rho1__ = rho10_ * g[1] + rho11_ * f[1];

    // Interpolate along axis 1 (x) → final density
    let rho_val = rho0__ * g[0] + rho1__ * f[0];

    // Gradient in lattice coordinates (partial derivatives)
    let grad_lat_0 = rho1__ - rho0__;  // ∂ρ/∂f1
    let grad_lat_1 = (rho10_ - rho00_) * g[0] + (rho11_ - rho01_) * f[0]; // ∂ρ/∂f2
    let rho_z0 = rho000 * g[1] * g[0] + rho010 * f[1] * g[0] + rho100 * g[1] * f[0] + rho110 * f[1] * f[0];
    let rho_z1 = rho001 * g[1] * g[0] + rho011 * f[1] * g[0] + rho101 * g[1] * f[0] + rho111 * f[1] * f[0];
    let grad_lat_2 = rho_z1 - rho_z0; // ∂ρ/∂f3

    // Transform to Cartesian gradient: grad_car = car2lat^T @ grad_lat
    let grad_lat = Vector3::new(grad_lat_0, grad_lat_1, grad_lat_2);

    (rho_val, chg.car2lat.transpose() * grad_lat)
}

/// Find the nearest grid point to continuous lattice coordinate `r`.
fn to_lat(chg: &ChargeGrid, r: [f64; 3]) -> [usize; 3] {
    // Round to nearest integer, then PBC fold
    let p = [
        r[0].round() as i32,
        r[1].round() as i32,
        r[2].round() as i32,
    ];
    chg.pbc_i(p)
}

/// One step of off-grid gradient ascent.  Returns `true` if moved.
fn step_offgrid(chg: &ChargeGrid, r: &mut [f64; 3], stepsize: f64) -> bool {
    let (_rho, grad_car) = rho_grad(chg, *r);
    let grad_norm = grad_car.norm();
    if grad_norm < 1e-30 {
        return false;
    }

    // Step in Cartesian, convert back to lattice
    let dr_car = grad_car * (stepsize / grad_norm);
    let dr_lat = chg.car2lat * dr_car;

    let r_new = [r[0] + dr_lat.x, r[1] + dr_lat.y, r[2] + dr_lat.z];
    let pm = to_lat(chg, r_new);
    let p_old = to_lat(chg, *r);

    // Check if density decreased at the new nearest grid point
    let rho_old = chg.rho_val([p_old[0] as i32, p_old[1] as i32, p_old[2] as i32]);
    let rho_new = chg.rho_val([pm[0] as i32, pm[1] as i32, pm[2] as i32]);

    if rho_new < rho_old {
        return false;
    }

    *r = r_new;
    true
}

/// Climb from continuous position `start` to a density maximum using off-grid steps.
fn max_offgrid(chg: &ChargeGrid, start: [f64; 3], stepsize: f64) -> [usize; 3] {
    let mut r = start;
    let max_iter = chg.nrho;
    for _ in 0..max_iter {
        if !step_offgrid(chg, &mut r, stepsize) {
            break;
        }
    }
    to_lat(chg, r)
}

/// Run Bader analysis using the off-grid method.
pub fn bader_offgrid(chg: &ChargeGrid, frame: &Frame, params: &BaderParams) -> BaderResult {
    let [n1, n2, n3] = chg.shape;
    let nrho = chg.nrho;
    let vol_cell = frame.cell.as_ref().unwrap().volume();

    // Default stepsize = minimum voxel dimension
    let stepsize = params.stepsize.unwrap_or_else(|| {
        let d0 = chg.lat2car.column(0).norm();
        let d1 = chg.lat2car.column(1).norm();
        let d2 = chg.lat2car.column(2).norm();
        d0.min(d1).min(d2)
    });

    let mut volnum = vec![0i32; nrho];
    let mut nvols = 0usize;
    let mut volpos_lat: Vec<[f64; 3]> = Vec::new();

    for i3 in 0..n3 {
        for i2 in 0..n2 {
            for i1 in 0..n1 {
                let idx = i1 + n1 * (i2 + n2 * i3);
                if volnum[idx] != 0 {
                    continue;
                }
                let start = [i1 as f64 + 0.5, i2 as f64 + 0.5, i3 as f64 + 0.5];
                let p_max = max_offgrid(chg, start, stepsize);
                let idx_max = p_max[0] + n1 * (p_max[1] + n2 * p_max[2]);

                let vol_id = if volnum[idx_max] > 0 {
                    volnum[idx_max]
                } else {
                    nvols += 1;
                    volpos_lat.push([p_max[0] as f64, p_max[1] as f64, p_max[2] as f64]);
                    nvols as i32
                };

                volnum[idx] = vol_id;
                if volnum[idx_max] == 0 {
                    volnum[idx_max] = vol_id;
                }
            }
        }
    }

    let _nvac = mark_vacuum(&mut volnum, &chg.rho, vol_cell, nrho, params.vacval);
    let mut volchg = calc_volchg(&volnum, &chg.rho, nvols, nrho);
    volchg.push(0.0);
    let (volpos_car, volpos_dir) = convert_volpos(&volpos_lat, chg, nvols);
    let cell = frame.cell.as_ref().unwrap();
    let (ionchg, iondist, nnion) = assign_chg2atom(&volpos_car, &volchg, frame, cell, nvols);
    let ionvol = calc_atomic_vol(&volnum, nvols, &nnion, frame.n_atoms(), vol_cell, nrho);
    let vacchg = if nvols < volchg.len() { volchg[nvols] } else { 0.0 };
    let vacvol = _nvac as f64 * vol_cell / nrho as f64;
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
