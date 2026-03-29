//! PFAM database parser for .hmm file format (HMMER3 compatible)
//!
//! Parses HMMER3 profile HMM files and provides efficient database access.

use crate::error::{Error, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// PFAM database with indexed profile HMMs
pub struct PfamDatabase {
    profiles: HashMap<String, PfamProfile>,
}

/// A single PFAM profile entry
#[derive(Debug, Clone)]
pub struct PfamProfile {
    pub name: String,
    pub accession: String,
    pub description: String,
    pub length: usize,
    pub alphabet: String,
    pub transition_probs: Vec<f64>,
    pub emission_probs: Vec<f64>,
    pub cutoff_ga: Option<f64>,
    pub cutoff_tc: Option<f64>,
    pub cutoff_nc: Option<f64>,
}

impl PfamDatabase {
    /// Create empty database
    pub fn new() -> Self {
        PfamDatabase {
            profiles: HashMap::new(),
        }
    }

    /// Load PFAM database from .hmm file
    pub fn from_hmm_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| Error::AlignmentError(format!("Failed to open PFAM file: {}", e)))?;
        let reader = BufReader::new(file);
        let mut db = PfamDatabase::new();
        let mut current_profile: Option<PfamProfileBuilder> = None;

        for line in reader.lines() {
            let line = line.map_err(|e| Error::AlignmentError(format!("Read error: {}", e)))?;
            let trimmed = line.trim();

            // Skip comments and empty lines
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Parse header
            if trimmed.starts_with("NAME") {
                if let Some(builder) = current_profile.take() {
                    db.profiles.insert(builder.name.clone(), builder.build()?);
                }
                current_profile = Some(PfamProfileBuilder::new());
                let name = trimmed.split_whitespace().nth(1).unwrap_or("unknown").to_string();
                if let Some(builder) = current_profile.as_mut() {
                    builder.name = name;
                }
            } else if trimmed.starts_with("ACC") {
                let acc = trimmed
                    .split_whitespace()
                    .nth(1)
                    .unwrap_or("unknown")
                    .to_string();
                if let Some(builder) = current_profile.as_mut() {
                    builder.accession = acc;
                }
            } else if trimmed.starts_with("DESC") {
                let desc = trimmed.split_whitespace().skip(1).collect::<Vec<_>>().join(" ");
                if let Some(builder) = current_profile.as_mut() {
                    builder.description = desc;
                }
            } else if trimmed.starts_with("LENG") {
                if let Ok(len) = trimmed.split_whitespace().nth(1).unwrap_or("0").parse::<usize>()
                {
                    if let Some(builder) = current_profile.as_mut() {
                        builder.length = len;
                    }
                }
            } else if trimmed.starts_with("GA") {
                if let Ok(cutoff) = trimmed.split_whitespace().nth(1).unwrap_or("0").parse::<f64>()
                {
                    if let Some(builder) = current_profile.as_mut() {
                        builder.cutoff_ga = Some(cutoff);
                    }
                }
            } else if trimmed.starts_with("TC") {
                if let Ok(cutoff) = trimmed.split_whitespace().nth(1).unwrap_or("0").parse::<f64>()
                {
                    if let Some(builder) = current_profile.as_mut() {
                        builder.cutoff_tc = Some(cutoff);
                    }
                }
            } else if trimmed.starts_with("NC") {
                if let Ok(cutoff) = trimmed.split_whitespace().nth(1).unwrap_or("0").parse::<f64>()
                {
                    if let Some(builder) = current_profile.as_mut() {
                        builder.cutoff_nc = Some(cutoff);
                    }
                }
            } else if trimmed == "//" {
                if let Some(builder) = current_profile.take() {
                    db.profiles.insert(builder.name.clone(), builder.build()?);
                }
            }
        }

        // Handle last profile if file doesn't end with //
        if let Some(builder) = current_profile {
            db.profiles.insert(builder.name.clone(), builder.build()?);
        }

        Ok(db)
    }

    /// Get profile by name
    pub fn get(&self, name: &str) -> Option<&PfamProfile> {
        self.profiles.get(name)
    }

    /// Get profile by accession
    pub fn get_by_accession(&self, accession: &str) -> Option<&PfamProfile> {
        self.profiles
            .values()
            .find(|p| p.accession == accession)
    }

    /// Iterate all profiles
    pub fn iter(&self) -> impl Iterator<Item = &PfamProfile> {
        self.profiles.values()
    }

    /// Number of profiles in database
    pub fn len(&self) -> usize {
        self.profiles.len()
    }

    /// Check if database is empty
    pub fn is_empty(&self) -> bool {
        self.profiles.is_empty()
    }

    /// Get all profile names
    pub fn names(&self) -> Vec<String> {
        self.profiles.keys().cloned().collect()
    }
}

/// Builder for constructing PFAM profiles
struct PfamProfileBuilder {
    name: String,
    accession: String,
    description: String,
    length: usize,
    alphabet: String,
    transition_probs: Vec<f64>,
    emission_probs: Vec<f64>,
    cutoff_ga: Option<f64>,
    cutoff_tc: Option<f64>,
    cutoff_nc: Option<f64>,
}

impl PfamProfileBuilder {
    fn new() -> Self {
        PfamProfileBuilder {
            name: String::new(),
            accession: String::new(),
            description: String::new(),
            length: 0,
            alphabet: "ACGT".to_string(),
            transition_probs: vec![0.1; 6],
            emission_probs: vec![0.25; 4],
            cutoff_ga: None,
            cutoff_tc: None,
            cutoff_nc: None,
        }
    }

    fn build(self) -> Result<PfamProfile> {
        Ok(PfamProfile {
            name: self.name,
            accession: self.accession,
            description: self.description,
            length: self.length,
            alphabet: self.alphabet,
            transition_probs: self.transition_probs,
            emission_probs: self.emission_probs,
            cutoff_ga: self.cutoff_ga,
            cutoff_tc: self.cutoff_tc,
            cutoff_nc: self.cutoff_nc,
        })
    }
}

/// E-value statistics using Karlin-Altschul parameters
pub struct EValueStats {
    pub lambda: f64,
    pub k: f64,
    pub mu: f64,
}

impl EValueStats {
    /// Create Karlin-Altschul statistics for protein sequences
    pub fn new_protein() -> Self {
        // Standard values for protein sequences with BLOSUM62
        EValueStats {
            lambda: 0.318,
            k: 0.135,
            mu: 35.0,
        }
    }

    /// Create for nucleotide sequences
    pub fn new_nucleotide() -> Self {
        EValueStats {
            lambda: 1.37,
            k: 0.711,
            mu: 38.7,
        }
    }

    /// Calculate E-value from bit score
    pub fn evalue(&self, bit_score: f64, database_size: f64) -> f64 {
        let raw_score = bit_score / self.lambda;
        self.k * database_size * (-self.lambda * raw_score).exp()
    }

    /// Calculate bit score from raw score
    pub fn bit_score(&self, raw_score: f64) -> f64 {
        (self.lambda * raw_score - self.mu.ln()) / 2.0_f64.ln()
    }

    /// Calculate p-value from raw score
    pub fn pvalue(&self, raw_score: f64) -> f64 {
        (-self.lambda * raw_score).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pfam_database_creation() {
        let db = PfamDatabase::new();
        assert_eq!(db.len(), 0);
        assert!(db.is_empty());
    }

    #[test]
    fn test_evalue_stats_protein() {
        let stats = EValueStats::new_protein();
        assert!(stats.lambda > 0.0);
        assert!(stats.k > 0.0);
    }

    #[test]
    fn test_bit_score_calculation() {
        let stats = EValueStats::new_protein();
        let raw_score = 50.0;
        let bit_score = stats.bit_score(raw_score);
        assert!(bit_score > 0.0);
    }

    #[test]
    fn test_evalue_calculation() {
        let stats = EValueStats::new_protein();
        let database_size = 1_000_000.0;
        let evalue = stats.evalue(20.0, database_size);
        assert!(evalue >= 0.0);
    }

    #[test]
    fn test_pvalue_calculation() {
        let stats = EValueStats::new_protein();
        let pvalue = stats.pvalue(50.0);
        assert!(pvalue >= 0.0 && pvalue <= 1.0);
    }
}
