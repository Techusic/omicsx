//! Batch sequence alignment with parallel processing
//!
//! Process multiple sequence alignments efficiently using thread pools.
//! Supports both local (Smith-Waterman) and global (Needleman-Wunsch) alignments.

use crate::alignment::{SmithWaterman, NeedlemanWunsch, AlignmentResult};
use crate::protein::Protein;
use crate::error::Result;
use rayon::prelude::*;

/// Configuration for batch alignment
#[derive(Debug, Clone)]
pub struct BatchConfig {
    /// Number of parallel threads (None = use Rayon default)
    pub num_threads: Option<usize>,
    
    /// Enable banded DP (with this bandwidth)
    pub bandwidth: Option<usize>,
    
    /// Use SIMD acceleration if available
    pub use_simd: bool,
}

impl Default for BatchConfig {
    fn default() -> Self {
        BatchConfig {
            num_threads: None,
            bandwidth: None,
            use_simd: false,
        }
    }
}

impl BatchConfig {
    /// Create a new batch configuration
    pub fn new() -> Self {
        Self::default()
    }

    /// Set number of threads
    pub fn with_threads(mut self, n: usize) -> Self {
        self.num_threads = Some(n);
        self
    }

    /// Enable banded DP
    pub fn with_bandwidth(mut self, bandwidth: usize) -> Self {
        self.bandwidth = Some(bandwidth);
        self
    }

    /// Enable SIMD
    pub fn with_simd(mut self, use_simd: bool) -> Self {
        self.use_simd = use_simd;
        self
    }
}

/// Item in a batch query - a single sequence to align against a reference
#[derive(Debug, Clone)]
pub struct BatchQuery {
    /// Query name/ID
    pub name: String,
    /// Sequence string
    pub sequence: String,
}

/// Result of a batch alignment
#[derive(Debug, Clone)]
pub struct BatchResult {
    /// Query name
    pub query_name: String,
    /// Alignment result
    pub alignment: AlignmentResult,
}

/// Batch Smith-Waterman local alignment
pub struct BatchSmithWaterman {
    config: BatchConfig,
    reference: Protein,
    aligner: SmithWaterman,
}

impl BatchSmithWaterman {
    /// Create a new batch Smith-Waterman aligner
    pub fn new(reference: &str, config: BatchConfig) -> Result<Self> {
        let reference = Protein::from_string(reference)?;
        let mut aligner = SmithWaterman::new();

        if let Some(bandwidth) = config.bandwidth {
            aligner = aligner.with_bandwidth(bandwidth);
        }

        if config.use_simd {
            aligner = aligner.with_simd(true);
        }

        Ok(BatchSmithWaterman {
            config,
            reference,
            aligner,
        })
    }

    /// Align multiple query sequences against the reference
    pub fn align_batch(&self, queries: Vec<BatchQuery>) -> Result<Vec<BatchResult>> {
        let reference = &self.reference;
        let aligner = &self.aligner;

        let results: Vec<std::result::Result<BatchResult, String>> = if let Some(num_threads) = self.config.num_threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .unwrap()
                .install(|| {
                    queries
                        .into_par_iter()
                        .map(|query| {
                            let query_protein = Protein::from_string(&query.sequence)
                                .map_err(|e| format!("Failed to parse query {}: {}", query.name, e))?;
                            let alignment = aligner.align(&query_protein, reference)
                                .map_err(|e| format!("Alignment failed for {}: {}", query.name, e))?;
                            Ok(BatchResult {
                                query_name: query.name,
                                alignment,
                            })
                        })
                        .collect()
                })
        } else {
            queries
                .into_par_iter()
                .map(|query| {
                    let query_protein = Protein::from_string(&query.sequence)
                        .map_err(|e| format!("Failed to parse query {}: {}", query.name, e))?;
                    let alignment = aligner.align(&query_protein, reference)
                        .map_err(|e| format!("Alignment failed for {}: {}", query.name, e))?;
                    Ok(BatchResult {
                        query_name: query.name,
                        alignment,
                    })
                })
                .collect()
        };

        // Convert Vec<Result> to Result<Vec>
        results.into_iter().collect::<std::result::Result<Vec<_>, String>>()
            .map_err(|e| crate::error::Error::Custom(e))
    }

    /// Filter alignment results by score threshold
    pub fn filter_by_score(results: &[BatchResult], min_score: i32) -> Vec<&BatchResult> {
        results
            .iter()
            .filter(|r| r.alignment.score >= min_score)
            .collect()
    }

    /// Filter alignment results by identity threshold
    pub fn filter_by_identity(results: &[BatchResult], min_identity: f64) -> Vec<&BatchResult> {
        results
            .iter()
            .filter(|r| r.alignment.identity() >= min_identity)
            .collect()
    }
}

/// Batch Needleman-Wunsch global alignment
pub struct BatchNeedlemanWunsch {
    config: BatchConfig,
    reference: Protein,
    aligner: NeedlemanWunsch,
}

impl BatchNeedlemanWunsch {
    /// Create a new batch Needleman-Wunsch aligner
    pub fn new(reference: &str, config: BatchConfig) -> Result<Self> {
        let reference = Protein::from_string(reference)?;
        let mut aligner = NeedlemanWunsch::new();

        if let Some(bandwidth) = config.bandwidth {
            aligner = aligner.with_bandwidth(bandwidth);
        }

        if config.use_simd {
            aligner = aligner.with_simd(true);
        }

        Ok(BatchNeedlemanWunsch {
            config,
            reference,
            aligner,
        })
    }

    /// Align multiple query sequences against the reference
    pub fn align_batch(&self, queries: Vec<BatchQuery>) -> Result<Vec<BatchResult>> {
        let reference = &self.reference;
        let aligner = &self.aligner;

        let results: Vec<std::result::Result<BatchResult, String>> = if let Some(num_threads) = self.config.num_threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .unwrap()
                .install(|| {
                    queries
                        .into_par_iter()
                        .map(|query| {
                            let query_protein = Protein::from_string(&query.sequence)
                                .map_err(|e| format!("Failed to parse query {}: {}", query.name, e))?;
                            let alignment = aligner.align(&query_protein, reference)
                                .map_err(|e| format!("Alignment failed for {}: {}", query.name, e))?;
                            Ok(BatchResult {
                                query_name: query.name,
                                alignment,
                            })
                        })
                        .collect()
                })
        } else {
            queries
                .into_par_iter()
                .map(|query| {
                    let query_protein = Protein::from_string(&query.sequence)
                        .map_err(|e| format!("Failed to parse query {}: {}", query.name, e))?;
                    let alignment = aligner.align(&query_protein, reference)
                        .map_err(|e| format!("Alignment failed for {}: {}", query.name, e))?;
                    Ok(BatchResult {
                        query_name: query.name,
                        alignment,
                    })
                })
                .collect()
        };

        // Convert Vec<Result> to Result<Vec>
        results.into_iter().collect::<std::result::Result<Vec<_>, String>>()
            .map_err(|e| crate::error::Error::Custom(e))
    }

    /// Filter alignment results by score threshold
    pub fn filter_by_score(results: &[BatchResult], min_score: i32) -> Vec<&BatchResult> {
        results
            .iter()
            .filter(|r| r.alignment.score >= min_score)
            .collect()
    }

    /// Filter alignment results by identity threshold
    pub fn filter_by_identity(results: &[BatchResult], min_identity: f64) -> Vec<&BatchResult> {
        results
            .iter()
            .filter(|r| r.alignment.identity() >= min_identity)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_batch_smith_waterman() -> Result<()> {
        let config = BatchConfig::new();
        let batch = BatchSmithWaterman::new("AGSGDSAF", config)?;

        let queries = vec![
            BatchQuery {
                name: "query1".to_string(),
                sequence: "AGSGD".to_string(),
            },
            BatchQuery {
                name: "query2".to_string(),
                sequence: "DSAF".to_string(),
            },
            BatchQuery {
                name: "query3".to_string(),
                sequence: "AGSGDSAF".to_string(),
            },
        ];

        let results = batch.align_batch(queries)?;
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].query_name, "query1");
        assert_eq!(results[1].query_name, "query2");
        assert_eq!(results[2].query_name, "query3");
        assert!(results[2].alignment.score > results[0].alignment.score); // Full match should score highest
        Ok(())
    }

    #[test]
    fn test_batch_needleman_wunsch() -> Result<()> {
        let config = BatchConfig::new();
        let batch = BatchNeedlemanWunsch::new("MGLSD", config)?;

        let queries = vec![
            BatchQuery {
                name: "q1".to_string(),
                sequence: "MGLS".to_string(),
            },
            BatchQuery {
                name: "q2".to_string(),
                sequence: "MGLSD".to_string(),
            },
        ];

        let results = batch.align_batch(queries)?;
        assert_eq!(results.len(), 2);
        assert!(results[1].alignment.score >= results[0].alignment.score);
        Ok(())
    }

    #[test]
    fn test_batch_with_bandwidth() -> Result<()> {
        let config = BatchConfig::new().with_bandwidth(10);
        let batch = BatchSmithWaterman::new("AGSGDSAFGCRESDVLQ", config)?;

        let queries = vec![
            BatchQuery {
                name: "similar".to_string(),
                sequence: "AGSGDSAFGCRESDVLQ".to_string(), // Identical
            },
        ];

        let results = batch.align_batch(queries)?;
        assert_eq!(results.len(), 1);
        assert!(results[0].alignment.score > 0);
        Ok(())
    }

    #[test]
    fn test_batch_filter_by_score() -> Result<()> {
        let config = BatchConfig::new();
        let batch = BatchSmithWaterman::new("AGSGDSAF", config)?;

        let queries = vec![
            BatchQuery {
                name: "high".to_string(),
                sequence: "AGSGDSAF".to_string(),
            },
            BatchQuery {
                name: "low".to_string(),
                sequence: "GGGGGGGG".to_string(),
            },
        ];

        let results = batch.align_batch(queries)?;
        let high_score = BatchSmithWaterman::filter_by_score(&results, 20);
        
        // The high-scoring alignment should be included
        assert!(high_score.iter().any(|r| r.query_name == "high"));
        Ok(())
    }
}
