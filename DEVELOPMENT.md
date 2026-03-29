# Development Guide

**Maintained by**: Raghav Maheshwari (@techusic)  
**Email**: raghavmkota@gmail.com  
**Repository**: https://github.com/techusic/omnics-x

This document provides detailed instructions for developing Omnics-X.

## Quick Start

```bash
# Clone and setup
git clone https://github.com/techusic/omnics-x.git
cd omnics-x
cargo build --release

# Run tests
cargo test --lib

# Run examples
cargo run --example basic_alignment --release
cargo run --example neon_alignment --release
cargo run --example bam_format --release
```

## Build Variants

```bash
# Debug build (faster compilation)
cargo build

# Release build (optimized)
cargo build --release

# Check without building
cargo check

# Clean build
cargo clean && cargo build --release
```

## Testing

```bash
# All tests
cargo test --lib

# Specific test
cargo test --lib protein::tests::test_amino_acid_from_code

# With output
cargo test --lib -- --nocapture

# Single-threaded (helps with debugging)
cargo test --lib -- --test-threads=1

# Ignored tests only
cargo test --lib -- --ignored

# Benchmarks
cargo bench --bench alignment_benchmarks -- --verbose
```

## Code Quality

```bash
# Format check
cargo fmt --check

# Format (apply fixes)
cargo fmt

# Lint with Clippy
cargo clippy --release

# Strict linting
cargo clippy --release -- -D warnings

# Documentation check
cargo doc --no-deps --document-private-items
```

## Performance Debugging

### Profile with Perf (Linux)

```bash
# Install perf if needed
sudo apt-get install linux-tools-generic

# Build with profiling info
cargo build --release

# Run with perf
perf record -g target/release/example_name
perf report
```

### Benchmark with Criterion

```bash
cargo bench --bench alignment_benchmarks -- --verbose --save-baseline main

# Compare against baseline
cargo bench -- --baseline main
```

### Check SIMD Instructions

```bash
# Disassemble to see SIMD instructions
objdump -d target/release/libomics_simd.so | grep -E "vmov|vpadd|vpmax"

# For macOS
otool -tV target/release/libomics_simd.dylib | grep -E "vmov|vpadd"
```

## Architecture Overview

### Phase 1: Protein Module (`src/protein/mod.rs`)

- **20 amino acids** (standard IUPAC codes)
- **Ambiguous codes** (B, Z, X, etc.)
- **Metadata support** (ID, description)
- **Serialization** via Serde

```rust
// Adding a new amino acid variant:
pub enum AminoAcid {
    // ...existing variants...
    NewAcid,
}

impl AminoAcid {
    pub fn new_from_code(code: char) -> Result<Self> {
        match code {
            // ...existing cases...
            'N' => Ok(AminoAcid::NewAcid),
            _ => Err(Error::InvalidCode(code)),
        }
    }
}
```

### Phase 2: Scoring Module (`src/scoring/mod.rs`)

- **Matrices**: BLOSUM62 (24×24), framework for PAM/GONNET
- **Penalties**: Affine gap model (-open, -extend)
- **Validation**: Dimension checks, value constraints

```rust
// Adding a new matrix:
impl ScoringMatrix {
    fn new_matrix_data() -> Vec<Vec<i32>> {
        vec![/* 24x24 matrix */]
    }
    
    pub fn new(matrix_type: MatrixType) -> Result<Self> {
        match matrix_type {
            MatrixType::NewMatrix => 
                Self::from_data(Self::new_matrix_data()),
            // ...existing cases...
        }
    }
}
```

### Phase 3: Alignment Module (`src/alignment/mod.rs`)

#### Kernels (`src/alignment/kernel/`)

- **scalar.rs**: Portable reference implementation
- **avx2.rs**: x86-64 8-wide parallelism
- **neon.rs**: ARM64 4-wide parallelism

```rust
// Kernel signature pattern:
fn smith_waterman_kernel(
    seq1_vec: &Vec<AminoAcid>,
    seq2_vec: &Vec<AminoAcid>,
    matrix: &ScoringMatrix,
    penalty: &AffinePenalty,
) -> AlignmentResult {
    // Implementation
}
```

#### Batch Processing (`src/alignment/batch.rs`)

- **Rayon parallelization** across queries
- **Filtering** by score/identity
- **Multi-threaded** execution

```rust
// Batch usage:
let batch = BatchSmithWaterman::new(reference, config)?;
let results = batch.align_batch(queries)?;
```

#### Binary Format (`src/alignment/bam.rs`)

- **Serialization** to binary
- **4-bit encoding** for sequences
- **CIGAR compression**

### Future Modules (`src/futures/`)

Each has clear `todo!()` markers for implementation:

1. **matrices.rs** - Additional scoring matrices
2. **formats.rs** - BLAST/GFF3 export
3. **gpu.rs** - GPU acceleration
4. **msa.rs** - Multiple sequence alignment
5. **hmm.rs** - Profile HMM
6. **phylogeny.rs** - Phylogenetic trees

## Testing Strategy

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feature() -> Result<()> {
        let input = /* setup */;
        let expected = /* known result */;
        assert_eq!(compute(input), expected);
        Ok(())
    }
}
```

### Integration Tests

Tests in `/tests/` directory test end-to-end workflows.

### Edge Cases to Test

- Empty sequences
- Single amino acid
- Very long sequences
- Sequences with ambiguous codes
- All same character
- Completely different sequences
- Matrix mismatches

## Commit Workflow

```bash
# Create feature branch
git checkout -b feature/my-feature

# Make changes
vim src/module/file.rs

# Stage changes
git add src/

# Commit with conventional format
git commit -m "feat(module): description"

# Push to fork
git push origin feature/my-feature

# Create PR on GitHub
```

## Debugging Tips

### Print Debugging

```rust
dbg!(variable);
eprintln!("Debug: {:?}", value);
```

### Run with Backtrace

```bash
RUST_BACKTRACE=1 cargo test --lib -- --nocapture
RUST_BACKTRACE=full cargo run --example basic_alignment
```

### Conditional Compilation

```rust
#[cfg(debug_assertions)]
eprintln!("Debug info");

#[cfg(all(test, feature = "debug"))]
fn dump_state() { /* ... */ }
```

## Performance Optimization Workflow

1. **Baseline**: `cargo bench --bench alignment_benchmarks -- --save-baseline main`
2. **Change code**: Make optimization
3. **Rerun**: `cargo bench --bench alignment_benchmarks -- --baseline main`
4. **Compare**: Review criterion output in `target/criterion/`

## Cross-Platform Testing

```bash
# Test on x86-64
cargo test --lib

# Build for ARM64 (if cross installed)
cargo build --release --target aarch64-unknown-linux-gnu

# Check SIMD features
rustc --print cfg | grep target_feature
```

## Documentation

```bash
# Generate and open docs
cargo doc --open

# Generate with private items
cargo doc --no-deps --document-private-items --open

# Check documentation compiles
cargo test --doc
```

## Common Issues

### "No such file or directory" errors
```bash
cargo clean
cargo build --release
```

### Test failures after git pull
```bash
git clean -fdx
cargo test --lib
```

### Benchmark time is unstable
- Close other applications
- Use criterion's statistical features
- Run multiple times with `--measurement-time`

### SIMD not vectorizing
```bash
rustc --print cfg | grep target_feature
cargo build --release --target x86_64-unknown-linux-gnu
```

## Release Workflow

See `.github/RELEASE_CHECKLIST.md` for detailed release process.

## Questions?

- Review existing tests for patterns
- Check `CONTRIBUTING.md` for guidelines
- See `README.md` for API examples
- Review `.github/copilot-instructions.md` for architecture notes
