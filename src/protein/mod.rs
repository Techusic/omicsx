//! # Protein Primitives Module
//!
//! Type-safe representations of amino acids and protein sequences.
//! This module provides the foundation for biological accuracy in sequence alignment.

use serde::{Deserialize, Serialize};
use std::fmt;
use crate::error::{Error, Result};

/// Standard amino acid types (20 canonical amino acids)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AminoAcid {
    Alanine,      // A
    Arginine,     // R
    Asparagine,   // N
    AsparticAcid, // D
    Cysteine,     // C
    GlutamicAcid, // E
    Glutamine,    // Q
    Glycine,      // G
    Histidine,    // H
    Isoleucine,   // I
    Leucine,      // L
    Lysine,       // K
    Methionine,   // M
    Phenylalanine,// F
    Proline,      // P
    Serine,       // S
    Threonine,    // T
    Tryptophan,   // W
    Tyrosine,     // Y
    Valine,       // V
    // Special codes
    Ambiguous,    // B (Aspartic or Asparagine)
    Uncertain,    // X (Unknown)
    Stop,         // * (Stop codon)
    Gap,          // - (Gap in alignment)
}

impl AminoAcid {
    /// Convert from single-letter IUPAC code to AminoAcid
    pub fn from_code(c: char) -> Result<Self> {
        match c.to_ascii_uppercase() {
            'A' => Ok(AminoAcid::Alanine),
            'R' => Ok(AminoAcid::Arginine),
            'N' => Ok(AminoAcid::Asparagine),
            'D' => Ok(AminoAcid::AsparticAcid),
            'C' => Ok(AminoAcid::Cysteine),
            'E' => Ok(AminoAcid::GlutamicAcid),
            'Q' => Ok(AminoAcid::Glutamine),
            'G' => Ok(AminoAcid::Glycine),
            'H' => Ok(AminoAcid::Histidine),
            'I' => Ok(AminoAcid::Isoleucine),
            'L' => Ok(AminoAcid::Leucine),
            'K' => Ok(AminoAcid::Lysine),
            'M' => Ok(AminoAcid::Methionine),
            'F' => Ok(AminoAcid::Phenylalanine),
            'P' => Ok(AminoAcid::Proline),
            'S' => Ok(AminoAcid::Serine),
            'T' => Ok(AminoAcid::Threonine),
            'W' => Ok(AminoAcid::Tryptophan),
            'Y' => Ok(AminoAcid::Tyrosine),
            'V' => Ok(AminoAcid::Valine),
            'B' => Ok(AminoAcid::Ambiguous),
            'X' => Ok(AminoAcid::Uncertain),
            '*' => Ok(AminoAcid::Stop),
            '-' => Ok(AminoAcid::Gap),
            _ => Err(Error::InvalidAminoAcid(c)),
        }
    }

    /// Get single-letter IUPAC code
    pub fn to_code(&self) -> char {
        match self {
            AminoAcid::Alanine => 'A',
            AminoAcid::Arginine => 'R',
            AminoAcid::Asparagine => 'N',
            AminoAcid::AsparticAcid => 'D',
            AminoAcid::Cysteine => 'C',
            AminoAcid::GlutamicAcid => 'E',
            AminoAcid::Glutamine => 'Q',
            AminoAcid::Glycine => 'G',
            AminoAcid::Histidine => 'H',
            AminoAcid::Isoleucine => 'I',
            AminoAcid::Leucine => 'L',
            AminoAcid::Lysine => 'K',
            AminoAcid::Methionine => 'M',
            AminoAcid::Phenylalanine => 'F',
            AminoAcid::Proline => 'P',
            AminoAcid::Serine => 'S',
            AminoAcid::Threonine => 'T',
            AminoAcid::Tryptophan => 'W',
            AminoAcid::Tyrosine => 'Y',
            AminoAcid::Valine => 'V',
            AminoAcid::Ambiguous => 'B',
            AminoAcid::Uncertain => 'X',
            AminoAcid::Stop => '*',
            AminoAcid::Gap => '-',
        }
    }

    /// Get unique index for matrix lookups (0-24 for standard + special codes)
    pub fn index(&self) -> usize {
        match self {
            AminoAcid::Alanine => 0,
            AminoAcid::Arginine => 1,
            AminoAcid::Asparagine => 2,
            AminoAcid::AsparticAcid => 3,
            AminoAcid::Cysteine => 4,
            AminoAcid::GlutamicAcid => 5,
            AminoAcid::Glutamine => 6,
            AminoAcid::Glycine => 7,
            AminoAcid::Histidine => 8,
            AminoAcid::Isoleucine => 9,
            AminoAcid::Leucine => 10,
            AminoAcid::Lysine => 11,
            AminoAcid::Methionine => 12,
            AminoAcid::Phenylalanine => 13,
            AminoAcid::Proline => 14,
            AminoAcid::Serine => 15,
            AminoAcid::Threonine => 16,
            AminoAcid::Tryptophan => 17,
            AminoAcid::Tyrosine => 18,
            AminoAcid::Valine => 19,
            AminoAcid::Ambiguous => 20,
            AminoAcid::Uncertain => 21,
            AminoAcid::Stop => 22,
            AminoAcid::Gap => 23,
        }
    }
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_code())
    }
}

/// A protein polymer (sequence of amino acids)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Protein {
    sequence: Vec<AminoAcid>,
    id: Option<String>,
    description: Option<String>,
}

impl Protein {
    /// Create a new protein from a sequence
    pub fn new(sequence: Vec<AminoAcid>) -> Result<Self> {
        if sequence.is_empty() {
            return Err(Error::EmptySequence);
        }
        Ok(Protein {
            sequence,
            id: None,
            description: None,
        })
    }

    /// Create a protein from IUPAC single-letter codes
    pub fn from_string(seq_str: &str) -> Result<Self> {
        let sequence: Result<Vec<_>> = seq_str
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(AminoAcid::from_code)
            .collect();
        
        Self::new(sequence?)
    }

    /// Set the protein ID
    pub fn with_id(mut self, id: String) -> Self {
        self.id = Some(id);
        self
    }

    /// Set the protein description
    pub fn with_description(mut self, description: String) -> Self {
        self.description = Some(description);
        self
    }

    /// Get the amino acid sequence
    pub fn sequence(&self) -> &[AminoAcid] {
        &self.sequence
    }

    /// Get the sequence length
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if sequence is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get protein ID
    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }

    /// Get protein description
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Get sequence as string representation
    pub fn as_string(&self) -> String {
        self.sequence.iter().map(|aa| aa.to_code()).collect()
    }
}

impl fmt::Display for Protein {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_amino_acid_from_code() {
        assert_eq!(AminoAcid::from_code('A').unwrap(), AminoAcid::Alanine);
        assert_eq!(AminoAcid::from_code('G').unwrap(), AminoAcid::Glycine);
        assert!(AminoAcid::from_code('Z').is_err());
    }

    #[test]
    fn test_amino_acid_to_code() {
        assert_eq!(AminoAcid::Alanine.to_code(), 'A');
        assert_eq!(AminoAcid::Glycine.to_code(), 'G');
    }

    #[test]
    fn test_protein_creation() {
        let seq = Protein::from_string("AGSGD").unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq.as_string(), "AGSGD");
    }

    #[test]
    fn test_protein_with_metadata() {
        let seq = Protein::from_string("AGS")
            .unwrap()
            .with_id("P123".to_string())
            .with_description("Test protein".to_string());
        
        assert_eq!(seq.id(), Some("P123"));
        assert_eq!(seq.description(), Some("Test protein"));
    }
}
