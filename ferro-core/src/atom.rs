//! Atom data: the owned [`Atom`] a frame is built from, the per-atom [`Site`] a frame
//! stores once and shares between frames, and the [`AtomRef`] view a frame hands out.

use nalgebra::Vector3;

/// One atom, owned: what a reader or a structure operation builds a frame from.
///
/// A frame does not store `Atom`s. It keeps what stays the same from frame to frame
/// (element, label, mass, magnetic moment: a [`Site`]) in one table shared between
/// frames, and what changes (position, charge) in arrays of its own; see
/// [`crate::Frame`]. An atom's identity is its index in the frame; no index field
/// is stored. When `mass` is `None`, [`Atom::effective_mass`] looks the element up.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Element symbol, e.g. "Fe", "O"
    pub element: String,
    /// Cartesian position \[Å\]
    pub position: Vector3<f64>,
    /// Optional site label, e.g. "Fe1", "Fe2", telling sites of one element apart
    pub label: Option<String>,
    /// Mass override \[amu\]; `None` means the element table's standard value
    pub mass: Option<f64>,
    /// Initial magnetic moment (the guess for a spin-polarised DFT run) \[μ_B\]
    pub magmom: Option<f64>,
    /// Atomic charge \[e\]: a LAMMPS `q` column, or a Bader / DDEC result
    pub charge: Option<f64>,
}

impl Atom {
    /// The minimal constructor: element symbol and position.
    pub fn new(element: impl Into<String>, position: Vector3<f64>) -> Self {
        Self {
            element: element.into(),
            position,
            label: None,
            mass: None,
            magmom: None,
            charge: None,
        }
    }

    /// Effective mass \[amu\]: `mass` if set, otherwise the element table's, 1.0 if unknown.
    pub fn effective_mass(&self) -> f64 {
        effective_mass(&self.element, self.mass)
    }

    /// This atom as the view a frame hands out.
    pub fn view(&self) -> AtomRef<'_> {
        AtomRef {
            element: &self.element,
            position: self.position,
            label: self.label.as_deref(),
            mass: self.mass,
            magmom: self.magmom,
            charge: self.charge,
        }
    }

    /// Cartesian distance to another atom \[Å\], ignoring periodicity.
    /// Periodic distances go through [`crate::periodic::PeriodicCell::minimum_image`].
    pub fn distance_to(&self, other: &Atom) -> f64 {
        (self.position - other.position).norm()
    }
}

fn effective_mass(element: &str, mass: Option<f64>) -> f64 {
    mass.unwrap_or_else(|| {
        crate::data::elements::by_symbol(element)
            .map(|e| e.atomic_mass)
            .unwrap_or(1.0)
    })
}

/// What an atom is, as opposed to where it is: the part a trajectory repeats in every
/// frame, so frames share one table of them (review B15).
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub element: String,
    pub label: Option<String>,
    pub mass: Option<f64>,
    pub magmom: Option<f64>,
}

impl Site {
    /// Effective mass \[amu\], as [`Atom::effective_mass`].
    pub fn effective_mass(&self) -> f64 {
        effective_mass(&self.element, self.mass)
    }
}

/// One atom of a frame, read out of the frame's site table and arrays.
///
/// Fields are copies or borrows, so `a.position - b.position` and
/// `a.element == "O"` read as they did on an owned [`Atom`].
#[derive(Debug, Clone, Copy)]
pub struct AtomRef<'a> {
    pub element: &'a str,
    pub position: Vector3<f64>,
    pub label: Option<&'a str>,
    pub mass: Option<f64>,
    pub magmom: Option<f64>,
    pub charge: Option<f64>,
}

impl AtomRef<'_> {
    /// Effective mass \[amu\], as [`Atom::effective_mass`].
    pub fn effective_mass(&self) -> f64 {
        effective_mass(self.element, self.mass)
    }

    /// An owned copy, for building another frame from this atom.
    pub fn to_atom(&self) -> Atom {
        Atom {
            element: self.element.to_string(),
            position: self.position,
            label: self.label.map(str::to_string),
            mass: self.mass,
            magmom: self.magmom,
            charge: self.charge,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let atom = Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(atom.element, "Fe");
        assert!(atom.label.is_none());
        assert!(atom.mass.is_none());
    }

    #[test]
    fn test_effective_mass_from_table() {
        let atom = Atom::new("Fe", Vector3::zeros());
        assert!((atom.effective_mass() - 55.845).abs() < 1e-3);
    }

    #[test]
    fn test_effective_mass_override() {
        let mut atom = Atom::new("Fe", Vector3::zeros());
        atom.mass = Some(56.0);
        assert!((atom.effective_mass() - 56.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance() {
        let a = Atom::new("O", Vector3::new(0.0, 0.0, 0.0));
        let b = Atom::new("H", Vector3::new(1.0, 0.0, 0.0));
        assert!((a.distance_to(&b) - 1.0).abs() < 1e-10);
    }
}
