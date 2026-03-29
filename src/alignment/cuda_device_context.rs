//! CUDA device context (legacy interface)
//!
//! This module is deprecated in favor of `cuda_runtime.rs` which provides
//! better cudarc integration. This file is kept for backward compatibility.

use crate::error::Result;

/// Deprecated: Use `GpuRuntime` from `cuda_runtime.rs` instead
#[derive(Debug, Clone)]
pub struct CudaDeviceContext {
    #[allow(dead_code)]
    device_id: u32,
}

impl CudaDeviceContext {
    /// Deprecated: Create a new CUDA device context
    /// Use `GpuRuntime::new()` instead
    pub fn new(device_id: u32) -> Result<Self> {
        Ok(CudaDeviceContext { device_id })
    }

    /// Deprecated: Get device ID
    pub fn device_id(&self) -> u32 {
        self.device_id
    }
}
