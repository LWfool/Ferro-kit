//! Atom data: the owned [`Atom`] a frame is built from, the per-atom [`Site`] a frame
//! stores once and shares between frames, and the [`AtomRef`] view a frame hands out.

use nalgebra::Vector3;

use crate::error::{ChemError, Result};

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

    /// Effective mass \[amu\]: `mass` if set, otherwise the element table's; `Err`
    /// naming the element when it is in neither.
    pub fn effective_mass(&self) -> Result<f64> {
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
}

/// No fallback (user decision 2026-10-20, review A-3): 1 amu used to stand in for an
/// unknown symbol, silently dragging densities down and skewing momenta, and guessing
/// from a prefix (`Fe1` -> Fe) would hide a typo. The user fixes the input instead.
fn effective_mass(element: &str, mass: Option<f64>) -> Result<f64> {
    if let Some(m) = mass {
        return Ok(m);
    }
    crate::data::elements::by_symbol(element).map(|e| e.atomic_mass).ok_or_else(|| {
        ChemError::ValidationError(format!(
            "no mass for element '{element}': not in the element table and no mass given; \
             correct the element name or give the mass"
        ))
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
    pub fn effective_mass(&self) -> Result<f64> {
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
    pub fn effective_mass(&self) -> Result<f64> {
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
        assert!((atom.effective_mass().unwrap() - 55.845).abs() < 1e-3);
    }

    #[test]
    fn test_effective_mass_override() {
        let mut atom = Atom::new("Fe", Vector3::zeros());
        atom.mass = Some(56.0);
        assert!((atom.effective_mass().unwrap() - 56.0).abs() < 1e-10);
    }

    // Review A-3: an unknown symbol used to count as 1 amu
    #[test]
    fn test_effective_mass_unknown_element_is_an_error() {
        let mut atom = Atom::new("Xx", Vector3::zeros());
        let e = atom.effective_mass().unwrap_err().to_string();
        assert!(e.contains("no mass for element 'Xx'"), "{e}");
        // No guessing from a prefix either
        assert!(Atom::new("Fe1", Vector3::zeros()).effective_mass().is_err());
        atom.mass = Some(3.0);
        assert_eq!(atom.effective_mass().unwrap(), 3.0);
    }
}
