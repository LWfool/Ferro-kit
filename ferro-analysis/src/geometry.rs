//! 几何分析工具：键长、键角、二面角、包围盒、回转半径

use ferro_core::Frame;
use nalgebra::Vector3;

/// 计算两原子间的 Cartesian 距离 \[Å\]（不考虑周期性）
pub fn distance(frame: &Frame, i: usize, j: usize) -> Option<f64> {
    let a = frame.atoms.get(i)?;
    let b = frame.atoms.get(j)?;
    Some((a.position - b.position).norm())
}

/// 计算三原子的键角 \[度\]，顶点为原子 j
pub fn angle(frame: &Frame, i: usize, j: usize, k: usize) -> Option<f64> {
    let ai = frame.atoms.get(i)?;
    let aj = frame.atoms.get(j)?;
    let ak = frame.atoms.get(k)?;
    let v1 = ai.position - aj.position;
    let v2 = ak.position - aj.position;
    let cos_a = (v1.dot(&v2) / (v1.norm() * v2.norm())).clamp(-1.0, 1.0);
    Some(cos_a.acos().to_degrees())
}

/// 计算四原子的二面角 \[度\]
pub fn dihedral(frame: &Frame, i: usize, j: usize, k: usize, l: usize) -> Option<f64> {
    let ai = frame.atoms.get(i)?;
    let aj = frame.atoms.get(j)?;
    let ak = frame.atoms.get(k)?;
    let al = frame.atoms.get(l)?;
    let b1 = aj.position - ai.position;
    let b2 = ak.position - aj.position;
    let b3 = al.position - ak.position;
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let m1 = n1.cross(&b2.normalize());
    Some(m1.dot(&n2).atan2(n1.dot(&n2)).to_degrees())
}

/// 计算回转半径 \[Å\]
pub fn radius_of_gyration(frame: &Frame) -> f64 {
    if frame.atoms.is_empty() { return 0.0; }
    let com = frame.center_of_mass();
    let total_mass = frame.total_mass();
    if total_mass == 0.0 { return 0.0; }
    let sum: f64 = frame
        .atoms
        .iter()
        .map(|a| {
            let r = (a.position - com).norm();
            a.effective_mass() * r * r
        })
        .sum();
    (sum / total_mass).sqrt()
}

/// 计算 Cartesian 包围盒尺寸 \[Å\]
pub fn bounding_box(frame: &Frame) -> Option<Vector3<f64>> {
    if frame.atoms.is_empty() { return None; }
    let first = frame.atoms[0].position;
    let (min, max) = frame.atoms.iter().skip(1).fold((first, first), |(mn, mx), a| {
        let p = a.position;
        (
            Vector3::new(mn.x.min(p.x), mn.y.min(p.y), mn.z.min(p.z)),
            Vector3::new(mx.x.max(p.x), mx.y.max(p.y), mx.z.max(p.z)),
        )
    });
    Some(max - min)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Atom;

    #[test]
    fn test_distance() {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("C", Vector3::new(1.0, 0.0, 0.0)));
        let d = distance(&frame, 0, 1).unwrap();
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_linear() {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O", Vector3::new(-1.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("O", Vector3::new(1.0, 0.0, 0.0)));
        let a = angle(&frame, 0, 1, 2).unwrap();
        assert!((a - 180.0).abs() < 1e-6);
    }

    #[test]
    fn test_bounding_box() {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("C", Vector3::new(2.0, 3.0, 1.0)));
        let bb = bounding_box(&frame).unwrap();
        assert!((bb.x - 2.0).abs() < 1e-10);
        assert!((bb.y - 3.0).abs() < 1e-10);
        assert!((bb.z - 1.0).abs() < 1e-10);
    }
}
