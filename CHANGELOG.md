# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-03-29

### Added - Phase 1: Protein Primitives
- Type-safe `AminoAcid` enum with IUPAC codes
- `Protein` struct with metadata support (ID, description, references)
- Serialization support via Serde (JSON, bincode)
- Comprehensive string parsing and validation

### Added - Phase 2: Scoring Infrastructure
- `ScoringMatrix` with BLOSUM62 data (24×24)
- Framework for PAM30/PAM70 matrices
- `AffinePenalty` with validation (enforces negative values)
- SAM format output with CIGAR string generation
- Penalty profiles: default, strict, liberal modes

### Added - Phase 3: SIMD Kernels
- Smith-Waterman and Needleman-Wunsch algorithms
- Scalar (portable) baseline implementations
- AVX2 kernel with 8-wide parallelism (x86-64)
- NEON kernel with 4-wide parallelism (ARM64)
- Runtime CPU feature detection
- Automatic kernel selection based on hardware
- CIGAR string generation for SAM/BAM compatibility

### Added - Advanced Features
- Banded DP algorithm (O(k·n) for similar sequences, ~10x speedup)
- Batch alignment API with Rayon parallelization
- BAM binary format (serialization/deserialization, ~4x compression vs SAM)
- Full SAM format support with header management

### Added - Testing & Documentation
- 32 comprehensive unit tests (100% pass rate)
- Criterion.rs benchmarks comparing SIMD vs scalar
- 4 production-ready examples
- Complete API documentation
- Cross-platform support (x86-64, ARM64, Windows/Linux/macOS)

### Added - Future Enhancement Scaffolding
- Module structure for 6 planned features
- Error types and data structures pre-defined
- 33 placeholder tests for future development

## [Unreleased]

### Planned - Scoring Matrices
- Additional PAM matrices (PAM40, PAM70, PAM120)
- GONNET statistical matrix
- HOXD matrix families
- Custom matrix loading from files

### Planned - BLAST Formats
- XML output (NCBI schema compatible)
- JSON serialization
- Tabular format (outfmt 6)
- GFF3 format output
- FASTA export

### Planned - GPU Acceleration
- CUDA backend (NVIDIA)
- HIP backend (AMD/ROCm)
- Vulkan compute shaders
- Multi-GPU load balancing
- Device detection and memory management

### Planned - Multiple Sequence Alignment
- Progressive MSA (ClustalW-like)
- Guide tree construction (UPGMA, neighbor-joining)
- Iterative refinement
- Profile alignment
- Consensus generation

### Planned - Profile HMM
- Hidden Markov model state machines
- Viterbi algorithm
- Forward-backward algorithm
- Domain detection
- Baum-Welch parameter optimization

### Planned - Phylogenetics
- UPGMA tree building
- Neighbor-joining algorithm
- Maximum parsimony
- Maximum likelihood inference
- Bootstrap resampling
- Newick format I/O

## Versioning

### 0.1.x
- Bug fixes and minor improvements
- No API breaking changes

### 0.2.0 (Next)
- Scoring Matrices enhancement
- BLAST format support

### 0.3.0
- GPU acceleration support

### 1.0.0
- All core features production-ready
- API stability guarantee
- Full Semantic Versioning commitment
