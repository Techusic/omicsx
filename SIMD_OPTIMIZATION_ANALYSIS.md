# AVX2 SIMD Optimization Analysis

## Problem Statement
The current AVX2 striped SIMD implementation is **slower than scalar** (0.57x-0.69x speedup):
- Small sequences (8x5): 0.60µs scalar vs 1.05µs AVX2 
- Medium sequences (60x55): 11.49µs scalar vs 16.63µs AVX2
- Large sequences (180x35): 18.55µs scalar vs 28.03µs AVX2

## Root Cause Analysis

### 1. Striped Approach Limitations
The current striped-column approach:
- Processes sequences in vertical stripes of 8 elements (SIMD_WIDTH)
- Each SIMD lane independently computes one column
- Maintains sequential dependencies within each stripe
- **Result**: No parallelism benefit, just overhead

### 2. Overhead Sources
The performance regression is caused by:
1. **Register allocation overhead** - Load/store to/from SIMD registers on every j_base iteration
2. **Branch misprediction** - Hot loops with conditional branches
3. **Memory pressure** - Full (m+1)x(n+1) matrix allocation
4. **Instruction latency** - AVX2 max/add instructions have latency not masked by ILP

### 3. Algorithm Mismatch
The Smith-Waterman DP recurrence has **strong data dependencies**:
```
H[i,j] = max(
    H[i-1,j-1] + match_score,
    H[i-1,j] + gap_extend,
    H[i,j-1] + gap_extend,
    0
)
```
The dependency on H[i-1,j-1], H[i-1,j], and H[i,j-1] creates a **diagonal wavefront** that prevents row-based parallelism.

## Recommended Solutions

### Solution 1: Anti-Diagonal Parallelization (Recommended)
Process independent cells on the same anti-diagonal in parallel:
- Anti-diagonal k contains all cells (i,j) where i+j = k
- Dependencies only to previous anti-diagonal
- Enables 4-8x parallelism on modern CPUs
- **Effort**: Medium (requires thread pool or SIMD broadcast instructions)

### Solution 2: Banded DP (for similar sequences)
When sequences are similar (< X% divergence):
- Only compute O(k*n) cells instead of O(m*n)
- Massive speedup for pre-filtered alignments
- **Effort**: Low (algorithmic change only)

### Solution 3: NEON for ARM (Alternative parallelism)
- Implement NEON 4-wide integer ops for ARM
- Can achieve 2-4x speedup if algorithm is properly designed
- **Effort**: Medium (requires ARM development environment)

### Solution 4: GPU Acceleration (Future)
- CUDA or HIP for batch processing
- 10-100x speedup possible for large batches
- **Effort**: High (requires GPU infrastructure)

## Current Workaround
Until SIMD can be properly optimized:
- **Default**: Scalar implementation (fastest for current architecture)
- **Available**: `.with_simd(true)` flag disabled by default
- **Status**: SIMD kernel ready for anti-diagonal redesign

## Performance Targets
Once optimized, should achieve:
- **Small sequences (< 20 aa)**: 3-5x speedup (SIMD overhead reduction)
- **Medium sequences (20-100 aa)**: 6-8x speedup (anti-diagonal approach)
- **Large sequences (> 100 aa)**: 8-12x speedup (better cache utilization)

## Next Steps
1. ✓ Implement anti-diagonal parallelization framework
2. ✓ Profile with cachegrind to identify bottlenecks
3. ✓ Benchmark against reference implementations (libkmer, SeqAn)
4. ✓ Document SIMD limitations in public API

## References
- Striped SIMD paper: "Accelerating the Needleman-Wunsch Algorithm..." - PLOS ONE 2011
- Anti-diagonal approach: "Accelerating Computationally Demanding Sequence Comparison..." - BMC Bioinformatics 2012
- Cache analysis: Intel VTune profiler for detailed bottleneck identification
