# OMICS-X: Comprehensive Technical Bug Audit
**Date**: April 1, 2026  
**Total Faults Identified**: 9 Critical/High Severity  
**Status**: Ready for Remediation

---

## Master Bug List

### ORIGINAL BUGS (4 Verified)

#### 1. ❌ NEON Excessive Allocations
**Severity**: 🟡 MEDIUM  
**File**: `src/alignment/simd_viterbi.rs` (lines 743-744)  
**Issue**: Heap allocation on every position causes O(N) allocation overhead  
**Impact**: ARM64 performance regression (Apple Silicon, AWS Graviton)  
**Fix**: Add scratch_neon_m and scratch_neon_backptr to ViterbiDecoder struct

#### 2. ❌ GPU Kernel Launch Placeholder (Viterbi)
**Severity**: 🔴 HIGH  
**File**: `src/alignment/simd_viterbi.rs` (line 260)  
**Issue**: `decode_cuda` has comment "would launch actual PTX kernel here" but runs CPU fallback  
**Impact**: GPU support advertised but computation happens on CPU  
**Fix**: Implement actual kernel launch with device.launch()

#### 3. ❌ GPU Kernel Launch Placeholder (smith_waterman/needleman_wunsch)
**Severity**: 🔴 HIGH  
**File**: `src/alignment/kernel/cuda.rs` (lines 166-180, 203+)  
**Issue**: CudaKernelManager returns CPU implementations despite GPU availability  
**Impact**: CUDA functions always use CPU code  
**Fix**: Implement actual GPU kernel execution

#### 4. ⚠️ HMMER3 Brittle Parsing
**Severity**: 🟡 MEDIUM  
**File**: `src/futures/hmmer3_full_parser.rs` (line 186)  
**Issue**: Uses `split_whitespace()` losing field alignment information  
**Impact**: Fragile to real-world HMMER3 files with special formatting  
**Fix**: Implement state-machine sub-parser for robust field extraction

---

### NEW BUGS (5 Verified + 1 Correct)

#### 5. ❌ Incomplete HMM Transition Parsing
**Severity**: 🔴 HIGH  
**File**: `src/alignment/hmmer3_parser.rs` (lines 410-453)  
**Issue**: Only reads 3 of 7 transitions (MM, MI, MD missing IM, II, DM, DD)  
**Impact**: Gaps cannot extend properly; multi-residue gaps biologically invalid  
**Fix**: Read all 7 transitions per state; update Hmmer3Model structure

#### 6. ❌ Hardcoded Shared Memory in CUDA
**Severity**: 🔴 **CRITICAL**  
**File**: `src/alignment/cuda_kernels_rtc.rs` (lines 131-132)  
**Issue**: Fixed `__shared__ float trans[512]` handles only 128 states max  
**Impact**: HMMs >1024 residues cause GPU page faults or corruption  
**Fix**: Use dynamic shared memory with `extern __shared__` and launch parameter

#### 7. ❌ Hardcoded Karlin-Altschul Statistics
**Severity**: 🔴 HIGH  
**File**: `src/futures/hmm.rs` (lines 570-572)  
**Issue**: Hardcoded λ=0.267, K=0.041; correct values λ=0.3176, K=0.134  
**Impact**: E-values off by 3-4×; researchers ignore real hits or trust false positives  
**Fix**: Parse LAMBDA/K from HMMER3 header or calculate from scoring matrix

#### 8. ❌ Missing I/O Compression Support
**Severity**: 🟡 MEDIUM  
**File**: `src/futures/cli_file_io.rs` (lines 75-88)  
**Issue**: No `.gz` file handling; plain `File::open` without GzDecoder  
**Impact**: Users must manually decompress multi-GB genomic databases  
**Fix**: Detect `.gz` extension and wrap reader with `flate2::GzDecoder`

#### 9. ❌ Unsafe Transition Indexing
**Severity**: 🔴 HIGH  
**File**: `src/alignment/simd_viterbi.rs` (lines 192-199) + `src/alignment/cuda_kernels_rtc.rs` (154-157)  
**Issue**: Assumes `transitions[0/1/2]` = MM/MI/MD without validation  
**Impact**: Silent corruption if transition order differs; combined with Bug #5, GPU gets wrong values  
**Fix**: Use explicit `TransitionType` enum mapping; dedicated `TransitionMatrix` struct

#### ✅ VERTICAL SIMD Optimization (NOT A BUG)
**Status**: VERIFIED CORRECT  
**File**: `src/alignment/simd_viterbi.rs` (lines 559-568)  
**Finding**: Already uses proper **Vertical/Striped SIMD** (processes states j, j+1, j+2, j+3 in parallel)  
**Action**: No fix needed

---

## Severity Breakdown

| Severity | Count | Issues |
|----------|-------|--------|
| 🔴 CRITICAL | 1 | Bug #6 (GPU memory overflow) |
| 🔴 HIGH | 5 | Bugs #2, #3, #5, #7, #9 |
| 🟡 MEDIUM | 3 | Bugs #1, #4, #8 |
| ✅ OK | 1 | SIMD optimization (correct) |

---

## Recommended Fix Order

1. **Bug #6** (Shared memory) - Highest impact, CRITICAL
2. **Bug #5** (Incomplete transitions) - Data corruption, affects all downstream code
3. **Bug #9** (Transition indexing) - Must be fixed after #5
4. **Bug #7** (Karlin-Altschul) - E-value correctness
5. **Bug #2, #3** (GPU kernels) - Implement actual GPU computation
6. **Bug #1** (NEON allocations) - Performance optimization
7. **Bug #4** (HMMER3 parsing) - Robustness improvement
8. **Bug #8** (Compression) - User experience

---

## Implementation Status

| Bug | Status | Estimated LOC |
|-----|--------|----------------|
| #1 | Pending | 15-20 |
| #2 | Pending | 40-50 |
| #3 | Pending | 80-100 |
| #4 | Pending | 50-70 |
| #5 | Pending | 60-80 |
| #6 | Pending | 30-40 |
| #7 | Pending | 20-30 |
| #8 | Pending | 40-50 |
| #9 | Pending | 70-90 |

**Total Estimated**: ~500-600 LOC changes

---

## Next Steps

- [ ] Fix Bug #6 (CUDA shared memory)
- [ ] Fix Bug #5 (HMM transitions)
- [ ] Fix Bug #9 (Transition indexing)
- [ ] Fix Bug #7 (Karlin-Altschul)
- [ ] Fix Bugs #2, #3 (GPU kernels)
- [ ] Fix Bug #1 (NEON allocations)
- [ ] Fix Bug #4 (HMMER3 parsing)
- [ ] Fix Bug #8 (Compression support)
- [ ] Run comprehensive test suite
- [ ] Deploy and validate

