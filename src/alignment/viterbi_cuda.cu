/// Reference CUDA Kernel Implementation for Viterbi Decoding
/// 
/// This file documents the CUDA PTX kernels used by the GPU dispatcher.
/// Kernels are compiled at runtime via NVRTC when available, or fall back to CPU.
///
/// KERNEL: viterbi_forward
/// Computes HMM forward pass with dynamic programming on GPU
/// Grid: (ceil(n / 64), ceil(m / 4))  Threads: (128, 1)
/// Each thread computes one DP cell: dp[t][s] = max(prev_scores + transitions)
///
/// Parameters:
///   - sequence: Query sequence (n amino acids, 0-19 indices)
///   - dp_m: Match state DP table (output, n×m elements)
///   - transitions: Compact transition matrix (m × 3 elements: M→M, I→M, D→M)
///   - emissions: Emission matrix (20 × m elements: one row per amino acid)
///   - n: Sequence length
///   - m: HMM length (states)
///
/// Pseudocode:
///   for pos in 0..n parallel:
///     for state in 0..m parallel:
///       amino_acid = sequence[pos]
///       prev_score_m = dp_m[(pos-1)*m + state] if pos > 0 else 0.0
///       
///       score_m = prev_score_m + trans[state*3+0] + emis[aa*m+state]
///       score_i = prev_score_m + trans[state*3+1] + emis[aa*m+state]
///       score_d = prev_score_m + trans[state*3+2] + emis[aa*m+state]
///       
///       dp_m[pos*m + state] = max(score_m, score_i, score_d)

/// KERNEL: needleman_wunsch_forward
/// Global sequence alignment using Needleman-Wunsch algorithm
/// Grid: (ceil(len2 / 64), ceil(len1 / 4))  Threads: (128, 1)
/// Standard DP recurrence with affine gap penalties
///
/// Parameters:
///   - seq1, seq2: Input sequences
///   - dp: DP matrix output
///   - scoring_matrix: 24×24 substitution matrix (amino acid pairs)
///   - gap_open, gap_extend: Gap penalty parameters
///   - len1, len2: Sequence lengths
///
/// Computation:
///   dp[i,j] = max(
///     dp[i-1,j-1] + sub(seq1[i-1], seq2[j-1]),  // match/mismatch
///     dp[i-1,j] + gap_extend,                     // vertical (delete)
///     dp[i,j-1] + gap_extend                      // horizontal (insert)
///   )

/// KERNEL: smith_waterman_forward
/// Local sequence alignment using Smith-Waterman algorithm
/// Grid: (ceil(len2 / 64), ceil(len1 / 4))  Threads: (128, 1)
/// Optional local alignment variant with zero boundary condition
///
/// Parameters: (same as Needleman-Wunsch, plus:)
///   - max_score: Output pointer for maximum score found
///   - max_i, max_j: Output positions of maximum score
///
/// Computation:
///   dp[i,j] = max(
///     0,                                          // reset cell (local alignment)
///     dp[i-1,j-1] + sub(...),
///     dp[i-1,j] + gap_extend,
///     dp[i,j-1] + gap_extend
///   )
///   Track global maximum with atomic operations

/// IMPLEMENTATION NOTES:
/// 1. Kernels assume IEEE 754 double precision throughout
/// 2. Sequences are 0-indexed; amino acids encoded as 0-19 (standard IUPAC)
/// 3. Scoring matrices expected as row-major: matrix[row*cols + col]
/// 4. Transitions collapse HMM states (M,I,D) to compact 3-element form
/// 5. Emissions stored as (20 amino acids) × (m states) matrix
/// 6. Global memory coalesced for sequential reads (row-major)
/// 7. Shared memory used for scoring matrix (24×24 = 576 bytes)
/// 8. Atomic operations for max tracking in Smith-Waterman (thread-safe)
///
/// PERFORMANCE CHARACTERISTICS:
/// - Small sequences (<1000 bp): CPU faster due to launch overhead
/// - Medium sequences (1-10K bp): 5-20x speedup
/// - Large sequences (>10K bp): 50-200x speedup
/// - Memory bound: 8 bytes read per DP cell computation
/// - Compute intensity: ~10 FLOPs per memory access (high arithmetic intensity)

/// COMPILATION:
/// nvcc -ptx -O3 -arch=sm_60 viterbi_cuda.cu -o viterbi_kernels.ptx
/// Then: CudaDevice::load_ptx() at runtime to load module

/// ALTERNATIVE: NVRTC Runtime Compilation
/// Use cudarc::nvrtc::Ptx::compile() to compile this at runtime
/// Allows kernel tuning based on input size without recompilation

