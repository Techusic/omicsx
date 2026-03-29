//! ARM NEON SIMD-accelerated alignment kernels for aarch64 architectures
//!
//! This module implements anti-diagonal parallelization of the Smith-Waterman and
//! Needleman-Wunsch algorithms using ARM NEON instructions.
//!
//! Key optimization: Process 4 parallel cells using NEON 32-bit integer vectors (i32x4).
//!
//! ## NEON Support
//!
//! ARM NEON is the Advanced SIMD extension for ARM processors:
//! - **Supported platforms**: ARMv7 (with NEON), ARMv8 (aarch64)
//! - **Vector width**: 128-bit (4x i32, 4x f32, 8x i16, 16x i8)
//! - **Common uses**: Mobile, embedded, server ARM CPUs (AWS Graviton, Apple Silicon)
//!
//! ## Implementation Strategy
//!
//! Similar to AVX2, uses anti-diagonal approach:
//! - Process 4 cells per NEON vector (i32x4)
//! - All cells on same anti-diagonal are independent
//! - Reduces complexity from striped column approach
//! - Target speedup: 2-4x over scalar for moderate-large sequences

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

use crate::protein::AminoAcid;
use crate::scoring::ScoringMatrix;
use crate::error::Result;

const NEON_WIDTH: usize = 4; // Number of i32 values in NEON register (128-bit / 32-bit)

/// NEON-optimized Smith-Waterman kernel using anti-diagonal approach
#[cfg(target_arch = "aarch64")]
pub fn smith_waterman_neon(
    seq1: &[AminoAcid],
    seq2: &[AminoAcid],
    matrix: &ScoringMatrix,
    _open_penalty: i32,
    extend_penalty: i32,
) -> Result<(Vec<Vec<i32>>, usize, usize)> {
    let m = seq1.len();
    let n = seq2.len();

    // Allocate DP matrix
    let mut h = vec![vec![0i32; n + 1]; m + 1];
    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    // Precompute scoring values
    let mut scores = vec![vec![0i32; n]; m];
    for i in 0..m {
        for j in 0..n {
            scores[i][j] = matrix.score(seq1[i], seq2[j]);
        }
    }

    // Safe wrappers for NEON intrinsics
    unsafe {
        let extend_vec = vdupq_n_s32(extend_penalty);
        let zero_vec = vdupq_n_s32(0);

        // Process anti-diagonals
        for k in 1..=(m + n) {
            let i_start = std::cmp::max(1, k as i32 - n as i32) as usize;
            let i_end = std::cmp::min(m, k - 1);

            if i_start > i_end {
                continue;
            }

            // Process cells in batches of NEON_WIDTH
            let mut batch_i = vec![0usize; NEON_WIDTH];
            let mut batch_j = vec![0usize; NEON_WIDTH];
            let mut diag_vals = vec![0i32; NEON_WIDTH];
            let mut up_vals = vec![0i32; NEON_WIDTH];
            let mut left_vals = vec![0i32; NEON_WIDTH];
            let mut scores_vals = vec![0i32; NEON_WIDTH];
            let mut results = vec![0i32; NEON_WIDTH];

            for batch_start in (i_start..=i_end).step_by(NEON_WIDTH) {
                let batch_end = std::cmp::min(batch_start + NEON_WIDTH, i_end + 1);
                let batch_len = batch_end - batch_start;

                // Prepare batch data
                for (batch_idx, i) in (batch_start..batch_end).enumerate() {
                    let j = k - i;
                    batch_i[batch_idx] = i;
                    batch_j[batch_idx] = j;

                    if i > 0 && j > 0 {
                        diag_vals[batch_idx] = h[i - 1][j - 1];
                    }
                    if i > 0 && j <= n {
                        up_vals[batch_idx] = h[i - 1][j];
                    }
                    if i <= m && j > 0 {
                        left_vals[batch_idx] = h[i][j - 1];
                    }
                    if i > 0 && j > 0 {
                        scores_vals[batch_idx] = scores[i - 1][j - 1];
                    }
                }

                // Pad unused slots
                for batch_idx in batch_len..NEON_WIDTH {
                    diag_vals[batch_idx] = 0;
                    up_vals[batch_idx] = 0;
                    left_vals[batch_idx] = 0;
                    scores_vals[batch_idx] = 0;
                }

                // Load into NEON registers
                let diag_vec = vld1q_s32(diag_vals.as_ptr());
                let up_vec = vld1q_s32(up_vals.as_ptr());
                let left_vec = vld1q_s32(left_vals.as_ptr());
                let scores_vec = vld1q_s32(scores_vals.as_ptr());

                // Compute branches
                let diag_result = vaddq_s32(diag_vec, scores_vec);
                let up_result = vaddq_s32(up_vec, extend_vec);
                let left_result = vaddq_s32(left_vec, extend_vec);

                // Max operations
                let max_du = vmaxq_s32(diag_result, up_result);
                let max_dul = vmaxq_s32(max_du, left_result);
                let result_vec = vmaxq_s32(max_dul, zero_vec);

                // Store results
                vst1q_s32(results.as_mut_ptr(), result_vec);

                // Update DP matrix and track maximum
                for (batch_idx, i) in (batch_start..batch_end).enumerate() {
                    let j = k - i;
                    h[i][j] = results[batch_idx];

                    if h[i][j] > max_score {
                        max_score = h[i][j];
                        max_i = i;
                        max_j = j;
                    }
                }
            }
        }
    }

    Ok((h, max_i, max_j))
}

/// Fallback for non-ARM architectures
#[cfg(not(target_arch = "aarch64"))]
pub fn smith_waterman_neon(
    _seq1: &[AminoAcid],
    _seq2: &[AminoAcid],
    _matrix: &ScoringMatrix,
    _open_penalty: i32,
    _extend_penalty: i32,
) -> Result<(Vec<Vec<i32>>, usize, usize)> {
    Err(crate::error::Error::Custom(
        "NEON kernel requires aarch64 architecture".to_string(),
    ))
}

/// NEON-optimized Needleman-Wunsch kernel
#[cfg(target_arch = "aarch64")]
pub fn needleman_wunsch_neon(
    seq1: &[AminoAcid],
    seq2: &[AminoAcid],
    matrix: &ScoringMatrix,
    open_penalty: i32,
    extend_penalty: i32,
) -> Result<Vec<Vec<i32>>> {
    let m = seq1.len();
    let n = seq2.len();

    let mut h = vec![vec![0i32; n + 1]; m + 1];

    // Initialize first row and column
    for i in 0..=m {
        h[i][0] = (i as i32) * open_penalty;
    }
    for j in 0..=n {
        h[0][j] = (j as i32) * open_penalty;
    }

    unsafe {
        let extend_vec = vdupq_n_s32(extend_penalty);

        // Process in batches similar to Smith-Waterman
        for i in 1..=m {
            for j_base in (1..=n).step_by(NEON_WIDTH) {
                let j_end = std::cmp::min(j_base + NEON_WIDTH, n + 1);
                let batch_len = j_end - j_base;

                let mut diag_vals = vec![0i32; NEON_WIDTH];
                let mut up_vals = vec![0i32; NEON_WIDTH];
                let mut left_vals = vec![0i32; NEON_WIDTH];
                let mut scores_vals = vec![0i32; NEON_WIDTH];
                let mut results = vec![0i32; NEON_WIDTH];

                for (idx, j) in (j_base..j_end).enumerate() {
                    diag_vals[idx] = h[i - 1][j - 1];
                    up_vals[idx] = h[i - 1][j];
                    left_vals[idx] = h[i][j - 1];
                    scores_vals[idx] = matrix.score(seq1[i - 1], seq2[j - 1]);
                }

                // Pad
                for idx in batch_len..NEON_WIDTH {
                    diag_vals[idx] = i32::MIN / 2;
                    up_vals[idx] = i32::MIN / 2;
                    left_vals[idx] = i32::MIN / 2;
                    scores_vals[idx] = 0;
                }

                // Load and compute
                let diag_vec = vld1q_s32(diag_vals.as_ptr());
                let up_vec = vld1q_s32(up_vals.as_ptr());
                let left_vec = vld1q_s32(left_vals.as_ptr());
                let scores_vec = vld1q_s32(scores_vals.as_ptr());

                let diag_result = vaddq_s32(diag_vec, scores_vec);
                let up_result = vaddq_s32(up_vec, extend_vec);
                let left_result = vaddq_s32(left_vec, extend_vec);

                let max_du = vmaxq_s32(diag_result, up_result);
                let result_vec = vmaxq_s32(max_du, left_result);

                vst1q_s32(results.as_mut_ptr(), result_vec);

                // Update matrix
                for (idx, j) in (j_base..j_end).enumerate() {
                    h[i][j] = results[idx];
                }
            }
        }
    }

    Ok(h)
}

/// Fallback for non-ARM architectures
#[cfg(not(target_arch = "aarch64"))]
pub fn needleman_wunsch_neon(
    _seq1: &[AminoAcid],
    _seq2: &[AminoAcid],
    _matrix: &ScoringMatrix,
    _open_penalty: i32,
    _extend_penalty: i32,
) -> Result<Vec<Vec<i32>>> {
    Err(crate::error::Error::Custom(
        "NEON kernel requires aarch64 architecture".to_string(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neon_smith_waterman_fallback() {
        let matrix = ScoringMatrix::default();
        let seq1 = vec![AminoAcid::Alanine, AminoAcid::Glycine];
        let seq2 = vec![AminoAcid::Alanine, AminoAcid::Glycine];

        let result = smith_waterman_neon(&seq1, &seq2, &matrix, -11, -1);
        if cfg!(target_arch = "aarch64") {
            assert!(result.is_ok(), "NEON should work on aarch64");
        } else {
            assert!(result.is_err(), "NEON should fail on non-ARM");
            assert_eq!(
                result.unwrap_err().to_string(),
                "NEON kernel requires aarch64 architecture"
            );
        }
    }

    #[test]
    fn test_neon_needleman_wunsch_fallback() {
        let matrix = ScoringMatrix::default();
        let seq1 = vec![AminoAcid::Alanine];
        let seq2 = vec![AminoAcid::Alanine];

        let result = needleman_wunsch_neon(&seq1, &seq2, &matrix, -11, -1);
        if cfg!(target_arch = "aarch64") {
            assert!(result.is_ok(), "NEON should work on aarch64");
        } else {
            assert!(result.is_err(), "NEON should fail on non-ARM");
            assert_eq!(
                result.unwrap_err().to_string(),
                "NEON kernel requires aarch64 architecture"
            );
        }
    }

    #[test]
    fn test_neon_kernel_marker() {
        // Marker test showing NEON kernel is available for aarch64 compilation
        // Real NEON acceleration requires actual aarch64 hardware (AWS Graviton, Apple Silicon, etc.)
        if cfg!(target_arch = "aarch64") {
            let matrix = ScoringMatrix::default();
            let seq1 = vec![
                AminoAcid::Alanine,
                AminoAcid::Glycine,
                AminoAcid::Valine,
                AminoAcid::Leucine,
            ];
            let seq2 = vec![
                AminoAcid::Alanine,
                AminoAcid::Glycine,
                AminoAcid::Valine,
                AminoAcid::Leucine,
            ];

            let (h, max_i, max_j) = smith_waterman_neon(&seq1, &seq2, &matrix, -11, -1)
                .expect("NEON should execute on aarch64");

            assert!(max_i > 0 && max_j > 0, "Should find maximum score position");
            assert!(h.len() > 0, "DP matrix should be populated");
            assert!(h.iter().all(|row| row.len() > 0), "All rows should be non-empty");
        }
    }
}
