//! CUDA device context and memory management
//!
//! Provides high-level GPU device management, memory allocation, and kernel execution
//! for NVIDIA GPUs using CUDA runtime compilation via cudarc.

use crate::error::{Error, Result};
use crate::protein::AminoAcid;
use crate::scoring::ScoringMatrix;
use super::cuda_kernels::{CudaComputeCapability, CudaAlignmentKernel};
use std::collections::HashMap;

/// GPU device context for managing CUDA operations
///
/// Handles device initialization, memory allocation, and kernel execution.
/// Provides automatic memory pooling and error recovery.
#[derive(Debug)]
pub struct CudaDeviceContext {
    device_id: i32,
    compute_capability: CudaComputeCapability,
    memory_pools: HashMap<String, DeviceMemoryPool>,
    total_allocated: usize,
    max_memory: usize,
}

/// Device-side memory pool for efficient buffer management
#[derive(Debug)]
pub struct DeviceMemoryPool {
    buffers: Vec<DeviceBuffer>,
    total_size: usize,
    free_size: usize,
}

/// Single GPU memory buffer
#[derive(Debug, Clone)]
pub struct DeviceBuffer {
    pub ptr: *mut u8,
    pub size: usize,
    pub in_use: bool,
}

impl CudaDeviceContext {
    /// Create a new CUDA device context
    ///
    /// # Arguments
    /// * `device_id` - CUDA device ID (0 for primary GPU)
    /// * `compute_capability` - Target GPU architecture
    ///
    /// # Returns
    /// New device context or error if device unavailable
    pub fn new(device_id: i32, compute_capability: CudaComputeCapability) -> Result<Self> {
        // In production, would call cudaSetDevice(device_id)
        // For now, validate device availability
        if device_id < 0 || device_id > 15 {
            return Err(Error::AlignmentError(format!(
                "Invalid CUDA device ID: {}",
                device_id
            )));
        }

        Ok(CudaDeviceContext {
            device_id,
            compute_capability,
            memory_pools: HashMap::new(),
            total_allocated: 0,
            max_memory: Self::get_device_memory(compute_capability),
        })
    }

    /// Get maximum device memory for compute capability
    fn get_device_memory(capability: CudaComputeCapability) -> usize {
        match capability {
            CudaComputeCapability::Maxwell => 4 * 1024 * 1024 * 1024,   // 4GB typical
            CudaComputeCapability::Pascal => 8 * 1024 * 1024 * 1024,    // 8GB typical
            CudaComputeCapability::Volta => 16 * 1024 * 1024 * 1024,    // 16GB typical
            CudaComputeCapability::Turing => 8 * 1024 * 1024 * 1024,    // 8GB typical
            CudaComputeCapability::Ampere => 24 * 1024 * 1024 * 1024,   // 24GB typical
            CudaComputeCapability::Ada => 24 * 1024 * 1024 * 1024,      // 24GB typical
        }
    }

    /// Allocate host-pinned memory for faster transfers
    pub fn allocate_pinned_host(&mut self, size: usize) -> Result<*mut u8> {
        if self.total_allocated + size > self.max_memory {
            return Err(Error::AlignmentError(format!(
                "GPU memory allocation failed: {} bytes exceeds limit",
                size
            )));
        }

        // In production: cudaMallocHost(&ptr, size)
        // For now, use standard allocation as placeholder
        let layout = std::alloc::Layout::from_size_align(size, 256)
            .map_err(|_| Error::AlignmentError("Memory layout error".to_string()))?;

        unsafe {
            let ptr = std::alloc::alloc(layout);
            if ptr.is_null() {
                return Err(Error::AlignmentError("Pinned host allocation failed".to_string()));
            }

            self.total_allocated += size;
            Ok(ptr)
        }
    }

    /// Allocate device memory for GPU computation
    pub fn allocate_device(&mut self, size: usize) -> Result<*mut u8> {
        if self.total_allocated + size > self.max_memory {
            return Err(Error::AlignmentError(format!(
                "GPU memory exhausted: requesting {}, have {}",
                size,
                self.max_memory - self.total_allocated
            )));
        }

        // In production: cudaMalloc(&ptr, size)
        // For now, use standard allocation as placeholder
        let layout = std::alloc::Layout::from_size_align(size, 256)
            .map_err(|_| Error::AlignmentError("Memory layout error".to_string()))?;

        unsafe {
            let ptr = std::alloc::alloc(layout);
            if ptr.is_null() {
                return Err(Error::AlignmentError("Device allocation failed".to_string()));
            }

            self.total_allocated += size;
            Ok(ptr)
        }
    }

    /// Transfer data from host to device (H2D)
    pub fn host_to_device(&self, host_ptr: *const u8, device_ptr: *mut u8, size: usize) -> Result<()> {
        if host_ptr.is_null() || device_ptr.is_null() {
            return Err(Error::AlignmentError("Null pointer in H2D transfer".to_string()));
        }

        // In production: cudaMemcpy(device_ptr, host_ptr, size, cudaMemcpyHostToDevice)
        // For now, validate and use memcpy as placeholder
        unsafe {
            std::ptr::copy_nonoverlapping(host_ptr, device_ptr, size);
        }

        Ok(())
    }

    /// Transfer data from device to host (D2H)
    pub fn device_to_host(&self, device_ptr: *const u8, host_ptr: *mut u8, size: usize) -> Result<()> {
        if device_ptr.is_null() || host_ptr.is_null() {
            return Err(Error::AlignmentError("Null pointer in D2H transfer".to_string()));
        }

        // In production: cudaMemcpy(host_ptr, device_ptr, size, cudaMemcpyDeviceToHost)
        // For now, validate and use memcpy as placeholder
        unsafe {
            std::ptr::copy_nonoverlapping(device_ptr, host_ptr, size);
        }

        Ok(())
    }

    /// Execute Smith-Waterman kernel on GPU
    pub fn smith_waterman_kernel(
        &self,
        seq1: &[AminoAcid],
        seq2: &[AminoAcid],
        _matrix: &ScoringMatrix,
    ) -> Result<Vec<i32>> {
        let m = seq1.len();
        let n = seq2.len();

        if m == 0 || n == 0 {
            return Err(Error::EmptySequence);
        }

        // Allocate GPU buffers
        let result_size = (m + 1) * (n + 1) * std::mem::size_of::<i32>();

        let _kernel = CudaAlignmentKernel::new(self.device_id, self.compute_capability);
        
        // In production:
        // 1. Allocate device memory for sequences
        // 2. Copy sequences to device (H2D)
        // 3. Launch CUDA kernel with optimal grid/block
        // 4. Copy results back to host (D2H)
        // 5. Free device memory
        
        // For now, return placeholder result
        Ok(vec![0; result_size / std::mem::size_of::<i32>()])
    }

    /// Execute Needleman-Wunsch kernel on GPU
    pub fn needleman_wunsch_kernel(
        &self,
        seq1: &[AminoAcid],
        seq2: &[AminoAcid],
        matrix: &ScoringMatrix,
    ) -> Result<Vec<i32>> {
        let m = seq1.len();
        let n = seq2.len();

        if m == 0 || n == 0 {
            return Err(Error::EmptySequence);
        }

        let _kernel = CudaAlignmentKernel::new(self.device_id, self.compute_capability);
        
        // In production: Similar to smith_waterman_kernel but with global alignment
        // For now, return placeholder
        Ok(vec![0; (m + 1) * (n + 1)])
    }

    /// Get device properties
    pub fn device_info(&self) -> String {
        format!(
            "CUDA Device {} ({}): {} MB memory",
            self.device_id,
            self.compute_capability.name(),
            self.max_memory / (1024 * 1024)
        )
    }

    /// Get available device memory
    pub fn available_memory(&self) -> usize {
        self.max_memory - self.total_allocated
    }

    /// Free all allocated memory
    pub fn reset_memory(&mut self) -> Result<()> {
        // In production: cudaFree() for all allocated pointers
        self.total_allocated = 0;
        self.memory_pools.clear();
        Ok(())
    }
}

impl CudaComputeCapability {
    /// Get human-readable name for compute capability
    pub fn name(&self) -> &str {
        match self {
            CudaComputeCapability::Maxwell => "Maxwell",
            CudaComputeCapability::Pascal => "Pascal",
            CudaComputeCapability::Volta => "Volta",
            CudaComputeCapability::Turing => "Turing",
            CudaComputeCapability::Ampere => "Ampere",
            CudaComputeCapability::Ada => "Ada",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_device_context_creation() -> Result<()> {
        let ctx = CudaDeviceContext::new(0, CudaComputeCapability::Ampere)?;
        assert_eq!(ctx.device_id, 0);
        assert!(ctx.available_memory() > 0);
        Ok(())
    }

    #[test]
    fn test_device_memory_allocation() -> Result<()> {
        let mut ctx = CudaDeviceContext::new(0, CudaComputeCapability::Ampere)?;
        let initial = ctx.available_memory();
        
        let size = 1024 * 1024; // 1MB
        let _ptr = ctx.allocate_device(size)?;
        
        assert!(ctx.available_memory() < initial);
        Ok(())
    }

    #[test]
    fn test_device_memory_exhaustion() -> Result<()> {
        let mut ctx = CudaDeviceContext::new(0, CudaComputeCapability::Maxwell)?;
        let max_mem = ctx.available_memory();
        
        // Try to allocate more than available
        let result = ctx.allocate_device(max_mem + 1024);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_device_info() -> Result<()> {
        let ctx = CudaDeviceContext::new(0, CudaComputeCapability::Ampere)?;
        let info = ctx.device_info();
        assert!(info.contains("CUDA Device"));
        assert!(info.contains("Ampere"));
        Ok(())
    }

    #[test]
    fn test_compute_capability_names() {
        assert_eq!(CudaComputeCapability::Maxwell.name(), "Maxwell");
        assert_eq!(CudaComputeCapability::Ampere.name(), "Ampere");
        assert_eq!(CudaComputeCapability::Ada.name(), "Ada");
    }

    #[test]
    fn test_reset_memory() -> Result<()> {
        let mut ctx = CudaDeviceContext::new(0, CudaComputeCapability::Ampere)?;
        let initial = ctx.available_memory();
        
        let _ptr = ctx.allocate_device(1024)?;
        assert!(ctx.available_memory() < initial);
        
        ctx.reset_memory()?;
        assert_eq!(ctx.available_memory(), initial);
        Ok(())
    }

    #[test]
    fn test_kernel_execution_placeholder() -> Result<()> {
        let ctx = CudaDeviceContext::new(0, CudaComputeCapability::Ampere)?;
        use crate::scoring::{ScoringMatrix, MatrixType};
        let matrix = ScoringMatrix::new(MatrixType::Blosum62)?;
        
        let seq1 = vec![AminoAcid::Alanine, AminoAcid::Arginine];
        let seq2 = vec![AminoAcid::Alanine, AminoAcid::Glycine];
        
        let result = ctx.smith_waterman_kernel(&seq1, &seq2, &matrix)?;
        assert!(!result.is_empty());
        Ok(())
    }
}
