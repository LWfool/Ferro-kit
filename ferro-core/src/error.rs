//! Error type definitions.

use thiserror::Error;

/// The error type of the library crates: the message alone. The category is named once, by
/// the layer that reports it (`FerroError`); a variant per cause used to stack a second and
/// third label on the same line, and nothing branched on it (review 2026-11-06 D8).
#[derive(Error, Debug)]
#[error("{0}")]
pub struct ChemError(pub String);

/// Result type alias
pub type Result<T> = std::result::Result<T, ChemError>;
