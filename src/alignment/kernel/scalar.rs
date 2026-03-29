//! Scalar (non-SIMD) alignment kernels for portability and validation

use crate::protein::AminoAcid;
use crate::scoring::ScoringMatrix;
use crate::error::Result;

/// Scalar Smith-Waterman kernel
/// 
/// Computes local alignment scores using standard dynamic programming.
/// This implementation serves as the baseline for SIMD optimization.
pub fn smith_waterman_scalar(
    seq1: &[AminoAcid],
    seq2: &[AminoAcid],
    matrix: &ScoringMatrix,
    _open_penalty: i32,
    extend_penalty: i32,
) -> Result<(Vec<Vec<i32>>, usize, usize)> {
    let m = seq1.len();
    let n = seq2.len();

    // Initialize DP matrix
    let mut h = vec![vec![0i32; n + 1]; m + 1];
    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    // Fill DP matrix with Smith-Waterman algorithm
    // Time complexity: O(m*n)
    for i in 1..=m {
        for j in 1..=n {
            let match_score = matrix.score(seq1[i - 1], seq2[j - 1]);
            let diagonal = h[i - 1][j - 1] + match_score;
            let up = h[i - 1][j] + extend_penalty;
            let left = h[i][j - 1] + extend_penalty;

            h[i][j] = std::cmp::max(0, std::cmp::max(diagonal, std::cmp::max(up, left)));

            if h[i][j] > max_score {
                max_score = h[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    Ok((h, max_i, max_j))
}

/// Scalar Needleman-Wunsch kernel
///
/// Computes global alignment scores using standard dynamic programming.
pub fn needleman_wunsch_scalar(
    seq1: &[AminoAcid],
    seq2: &[AminoAcid],
    matrix: &ScoringMatrix,
    open_penalty: i32,
    extend_penalty: i32,
) -> Result<Vec<Vec<i32>>> {
    let m = seq1.len();
    let n = seq2.len();

    // Initialize DP matrix with gap penalties for first row/column
    let mut h = vec![vec![0i32; n + 1]; m + 1];

    for i in 0..=m {
        h[i][0] = (i as i32) * open_penalty;
    }
    for j in 0..=n {
        h[0][j] = (j as i32) * open_penalty;
    }

    // Fill DP matrix with Needleman-Wunsch algorithm
    for i in 1..=m {
        for j in 1..=n {
            let match_score = matrix.score(seq1[i - 1], seq2[j - 1]);
            let diagonal = h[i - 1][j - 1] + match_score;
            let up = h[i - 1][j] + extend_penalty;
            let left = h[i][j - 1] + extend_penalty;

            h[i][j] = std::cmp::max(diagonal, std::cmp::max(up, left));
        }
    }

    Ok(h)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smith_waterman_scalar() -> Result<()> {
        let matrix = ScoringMatrix::default();
        let seq1 = vec![AminoAcid::Alanine, AminoAcid::Glycine, AminoAcid::Serine];
        let seq2 = vec![AminoAcid::Alanine, AminoAcid::Serine];

        let (h, max_i, max_j) = smith_waterman_scalar(&seq1, &seq2, &matrix, -11, -1)?;
        
        assert!(max_i > 0);
        assert!(max_j > 0);
        assert!(h[max_i][max_j] >= 0);
        
        Ok(())
    }
}
