use thiserror::Error;

/// Result type for omics-simd operations
pub type Result<T> = std::result::Result<T, Error>;

/// Errors that can occur in omics-simd
#[derive(Debug, Error)]
pub enum Error {
    #[error("Invalid amino acid code: {0}")]
    InvalidAminoAcid(char),

    #[error("Empty sequence")]
    EmptySequence,

    #[error("Sequence length mismatch")]
    LengthMismatch,

    #[error("Invalid scoring matrix dimensions")]
    InvalidMatrixDimensions,

    #[error("Invalid gap penalty")]
    InvalidGapPenalty,

    #[error("Alignment error: {0}")]
    AlignmentError(String),

    #[error("Serialization error: {0}")]
    SerializationError(#[from] bincode::Error),

    #[error("Custom error: {0}")]
    Custom(String),
}

impl From<String> for Error {
    fn from(msg: String) -> Self {
        Error::Custom(msg)
    }
}

impl From<&str> for Error {
    fn from(msg: &str) -> Self {
        Error::Custom(msg.to_string())
    }
}
