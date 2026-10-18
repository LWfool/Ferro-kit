//! Error type definitions.

use thiserror::Error;

/// The error type of the library crates.
#[derive(Error, Debug)]
pub enum ChemError {
    /// Text that does not parse
    #[error("Parse error: {0}")]
    ParseError(String),

    /// Data or parameters that do not hold
    #[error("Validation error: {0}")]
    ValidationError(String),

    /// An I/O error
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Result type alias
pub type Result<T> = std::result::Result<T, ChemError>;
