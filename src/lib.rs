//! # OMICS-SIMD: High-Performance Sequence Alignment Library
//!
//! A Rust library providing SIMD-accelerated sequence alignment algorithms
//! for petabyte-scale genomic data processing.
//!
//! ## Features
//!
//! - **Protein Primitives**: Type-safe amino acid and protein polymer representations
//! - **Scoring Infrastructure**: BLOSUM and PAM matrices with affine gap penalties
//! - **SIMD Kernels**: AVX2/NEON-optimized Smith-Waterman and Needleman-Wunsch implementations
//! - **Comprehensive Benchmarking**: Criterion.rs benchmarks comparing SIMD vs. scalar performance
//!
//! ## Project Phases
//!
//! ### Phase 1: Foundation
//! Extend the omics crate with amino acid types and protein polymers,
//! ensuring memory-safety and type-safety standards.
//!
//! ### Phase 2: Infrastructure
//! Develop scoring engine with affine gap penalties and CIGAR string generation.
//!
//! ### Phase 3: Optimization
//! Implement SIMD kernels using std::arch (AVX2, NEON) for high-performance computing.
//!
//! ## Example
//!
//! ```ignore
//! use omics_simd::protein::AminoAcid;
//! use omics_simd::alignment::SmithWaterman;
//!
//! // Create protein sequences
//! let seq1 = vec![AminoAcid::Alanine, AminoAcid::Glycine, AminoAcid::Serine];
//! let seq2 = vec![AminoAcid::Alanine, AminoAcid::Aspartic, AminoAcid::Serine];
//!
//! // Perform SIMD-accelerated alignment
//! let aligner = SmithWaterman::new();
//! let result = aligner.align(&seq1, &seq2);
//! ```

pub mod error;
pub mod protein;
pub mod scoring;
pub mod alignment;
pub mod futures;

pub use error::{Error, Result};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
