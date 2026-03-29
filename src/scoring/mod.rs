//! # Scoring Infrastructure Module
//!
//! Provides scoring matrices (BLOSUM, PAM) and affine gap penalty models
//! for biologically accurate sequence alignment.

use serde::{Deserialize, Serialize};
use crate::error::{Error, Result};
use crate::protein::AminoAcid;

/// Affine gap penalty model
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AffinePenalty {
    /// Gap opening penalty (typically negative)
    pub open: i32,
    /// Gap extension penalty (typically negative, less severe than open)
    pub extend: i32,
}

impl AffinePenalty {
    /// Create new affine penalty with validation
    pub fn new(open: i32, extend: i32) -> Result<Self> {
        if open > 0 || extend > 0 {
            return Err(Error::InvalidGapPenalty);
        }
        Ok(AffinePenalty { open, extend })
    }

    /// Default penalties suitable for protein alignment
    pub fn default_protein() -> Self {
        AffinePenalty {
            open: -11,
            extend: -1,
        }
    }

    /// Strict penalties for high-confidence alignment
    pub fn strict() -> Self {
        AffinePenalty {
            open: -16,
            extend: -4,
        }
    }

    /// Liberal penalties for distant sequences
    pub fn liberal() -> Self {
        AffinePenalty {
            open: -8,
            extend: -1,
        }
    }
}

impl Default for AffinePenalty {
    fn default() -> Self {
        Self::default_protein()
    }
}

/// Scoring matrix types available
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MatrixType {
    /// BLOSUM62 - most commonly used for protein alignment
    Blosum62,
    /// BLOSUM45 - for more distant sequences
    Blosum45,
    /// BLOSUM80 - for closely related sequences
    Blosum80,
    /// PAM30 - for distant sequences
    Pam30,
    /// PAM70 - for moderate divergence
    Pam70,
}

impl fmt::Display for MatrixType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixType::Blosum62 => write!(f, "BLOSUM62"),
            MatrixType::Blosum45 => write!(f, "BLOSUM45"),
            MatrixType::Blosum80 => write!(f, "BLOSUM80"),
            MatrixType::Pam30 => write!(f, "PAM30"),
            MatrixType::Pam70 => write!(f, "PAM70"),
        }
    }
}

/// Scoring matrix for amino acid substitutions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoringMatrix {
    matrix_type: MatrixType,
    scores: Vec<Vec<i32>>,
    size: usize,
}

impl ScoringMatrix {
    /// Create a new scoring matrix
    pub fn new(matrix_type: MatrixType) -> Result<Self> {
        let scores = match matrix_type {
            MatrixType::Blosum62 => Self::blosum62_data(),
            MatrixType::Blosum45 => Self::blosum45_data(),
            MatrixType::Blosum80 => Self::blosum80_data(),
            MatrixType::Pam30 => Self::pam30_data(),
            MatrixType::Pam70 => Self::pam70_data(),
        };

        let size = scores.len();
        if size != 24 {
            return Err(Error::InvalidMatrixDimensions);
        }

        Ok(ScoringMatrix {
            matrix_type,
            scores,
            size,
        })
    }

    /// Get score for amino acid pair
    pub fn score(&self, aa1: AminoAcid, aa2: AminoAcid) -> i32 {
        let i = aa1.index();
        let j = aa2.index();
        if i < self.size && j < self.size {
            self.scores[i][j]
        } else {
            -100 // Penalty for invalid indices
        }
    }

    /// Get matrix type
    pub fn matrix_type(&self) -> MatrixType {
        self.matrix_type
    }

    /// BLOSUM62 scoring matrix (most commonly used)
    fn blosum62_data() -> Vec<Vec<i32>> {
        vec![
            vec![4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -4], // A
            vec![-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4], // R
            vec![-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4], // N
            vec![-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4], // D
            vec![0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4], // C
            vec![-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 4, -1, -4], // E
            vec![-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4], // Q
            vec![0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -4, -3, -3, -1, -2, -1, -4], // G
            vec![-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4], // H
            vec![-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4], // I
            vec![-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4], // L
            vec![-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4], // K
            vec![-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4], // M
            vec![-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4], // F
            vec![-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4], // P
            vec![1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4], // S
            vec![0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4], // T
            vec![-3, -3, -4, -4, -2, -2, -3, -4, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4], // W
            vec![-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4], // Y
            vec![0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4], // V
            vec![-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4], // B
            vec![-1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 5, -1, -4], // Z
            vec![-1, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4], // X
            vec![-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1], // *
        ]
    }

    /// BLOSUM45 scoring matrix (for more distant sequences)
    fn blosum45_data() -> Vec<Vec<i32>> {
        // Placeholder implementation - use actual BLOSUM45 values
        Self::blosum62_data() // TODO: Replace with actual BLOSUM45
    }

    /// BLOSUM80 scoring matrix (for closely related sequences)
    fn blosum80_data() -> Vec<Vec<i32>> {
        // Placeholder implementation - use actual BLOSUM80 values
        Self::blosum62_data() // TODO: Replace with actual BLOSUM80
    }

    /// PAM30 scoring matrix
    fn pam30_data() -> Vec<Vec<i32>> {
        // Placeholder implementation - use actual PAM30 values
        Self::blosum62_data() // TODO: Replace with actual PAM30
    }

    /// PAM70 scoring matrix
    fn pam70_data() -> Vec<Vec<i32>> {
        // Placeholder implementation - use actual PAM70 values
        Self::blosum62_data() // TODO: Replace with actual PAM70
    }
}

impl Default for ScoringMatrix {
    fn default() -> Self {
        Self::new(MatrixType::Blosum62).expect("BLOSUM62 matrix should be valid")
    }
}

use std::fmt;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_affine_penalty() {
        let penalty = AffinePenalty::new(-11, -1).unwrap();
        assert_eq!(penalty.open, -11);
        assert_eq!(penalty.extend, -1);
    }

    #[test]
    fn test_invalid_penalty() {
        assert!(AffinePenalty::new(11, -1).is_err());
        assert!(AffinePenalty::new(-11, 1).is_err());
    }

    #[test]
    fn test_scoring_matrix() {
        let matrix = ScoringMatrix::default();
        let aa1 = AminoAcid::Alanine;
        let aa2 = AminoAcid::Alanine;
        assert_eq!(matrix.score(aa1, aa2), 4); // Diagonal should be positive
    }
}
