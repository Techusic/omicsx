//! # SIMD Kernel Module
//!
//! Architecture-specific SIMD-accelerated alignment kernels.
//! Provides compile-time selection between scalar and SIMD implementations.

pub mod scalar;
pub mod banded;

#[cfg(target_arch = "x86_64")]
pub mod avx2;

#[cfg(target_arch = "aarch64")]
pub mod neon;
