//! GPU runtime management with cudarc integration
//!
//! Provides safe wrappers around cudarc for memory management, kernel execution,
//! and device synchronization.

use crate::error::{Error, Result};
use std::sync::{Arc, Mutex};

#[cfg(feature = "cuda")]
use cudarc::driver::{CudaDevice, DevicePtr};

/// GPU device runtime with memory and kernel management
#[derive(Clone)]
pub struct GpuRuntime {
    #[cfg(feature = "cuda")]
    device: Arc<CudaDevice>,
    device_id: u32,
    total_memory: u64,
    allocated: Arc<Mutex<u64>>,
}

impl GpuRuntime {
    /// Initialize GPU runtime for all available devices
    pub fn detect_available_devices() -> Result<Vec<u32>> {
        #[cfg(feature = "cuda")]
        {
            use cudarc::driver::DriverError;
            match CudaDevice::list_devices() {
                Ok(devices) => Ok(devices.into_iter().map(|d| d as u32).collect()),
                Err(DriverError::NoCudaSupport | DriverError::NoDevices) => Ok(Vec::new()),
                Err(e) => Err(Error::AlignmentError(format!("CUDA device detection failed: {}", e))),
            }
        }

        #[cfg(not(feature = "cuda"))]
        {
            Ok(Vec::new())
        }
    }

    /// Create new GPU runtime for specified device
    pub fn new(device_id: u32) -> Result<Self> {
        #[cfg(feature = "cuda")]
        {
            use cudarc::driver::{CudaDevice, DriverError};
            let device = CudaDevice::new(device_id as usize)
                .map_err(|e| Error::AlignmentError(format!("Failed to initialize GPU {}: {}", device_id, e)))?;

            // Get device memory
            let total_memory = match device.get_device_memory() {
                Ok((free, total)) => total as u64,
                Err(e) => {
                    return Err(Error::AlignmentError(format!(
                        "Failed to query device memory: {}",
                        e
                    )))
                }
            };

            Ok(GpuRuntime {
                device: Arc::new(device),
                device_id,
                total_memory,
                allocated: Arc::new(Mutex::new(0)),
            })
        }

        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::AlignmentError(
                "CUDA support not compiled (enable 'cuda' feature)".to_string(),
            ))
        }
    }

    /// Get device ID
    pub fn device_id(&self) -> u32 {
        self.device_id
    }

    /// Get total device memory in bytes
    pub fn total_memory(&self) -> u64 {
        self.total_memory
    }

    /// Get allocated memory in bytes
    pub fn allocated_memory(&self) -> u64 {
        *self.allocated.lock().unwrap()
    }

    /// Get available memory in bytes
    pub fn available_memory(&self) -> u64 {
        self.total_memory - self.allocated_memory()
    }

    /// Allocate device memory
    pub fn allocate<T>(&self, size: usize) -> Result<GpuBuffer<T>>
    where
        T: Default + Clone + std::marker::Send,
    {
        let _byte_size = size * std::mem::size_of::<T>();

        #[cfg(feature = "cuda")]
        {
            let ptr = self.device
                .alloc::<T>(size)
                .map_err(|e| Error::AlignmentError(format!("GPU allocation failed: {}", e)))?;

            let mut allocated = self.allocated.lock().unwrap();
            *allocated += byte_size as u64;

            Ok(GpuBuffer {
                ptr: Some(ptr),
                size,
                device: Arc::clone(&self.device),
                allocated: Arc::clone(&self.allocated),
            })
        }

        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::AlignmentError(
                "CUDA support not available".to_string(),
            ))
        }
    }

    /// Copy data to device
    pub fn copy_to_device<T>(&self, _data: &[T]) -> Result<GpuBuffer<T>>
    where
        T: Default + Clone + std::marker::Send,
    {
        #[cfg(feature = "cuda")]
        {
            let buf = self.allocate::<T>(data.len())?;
            self.device
                .HTD::<T>(buf.ptr.as_ref().unwrap(), data)
                .map_err(|e| Error::AlignmentError(format!("H2D transfer failed: {}", e)))?;
            Ok(buf)
        }

        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::AlignmentError(
                "CUDA support not available".to_string(),
            ))
        }
    }

    /// Copy data from device to host
    pub fn copy_from_device<T>(&self, _buf: &GpuBuffer<T>) -> Result<Vec<T>>
    where
        T: Default + Clone + std::marker::Send,
    {
        #[cfg(feature = "cuda")]
        {
            let mut host_data = vec![T::default(); buf.size];
            self.device
                .DTH::<T>(host_data.as_mut_slice(), buf.ptr.as_ref().unwrap())
                .map_err(|e| Error::AlignmentError(format!("D2H transfer failed: {}", e)))?;
            Ok(host_data)
        }

        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::AlignmentError(
                "CUDA support not available".to_string(),
            ))
        }
    }

    /// Get device properties
    pub fn device_properties(&self) -> String {
        format!(
            "GPU {}: {} MB total, {} MB available",
            self.device_id,
            self.total_memory / (1024 * 1024),
            self.available_memory() / (1024 * 1024)
        )
    }

    /// Synchronize device (wait for all operations to complete)
    pub fn synchronize(&self) -> Result<()> {
        #[cfg(feature = "cuda")]
        {
            self.device
                .synchronize()
                .map_err(|e| Error::AlignmentError(format!("Device synchronization failed: {}", e)))
        }

        #[cfg(not(feature = "cuda"))]
        {
            Ok(())
        }
    }
}

/// GPU memory buffer wrapper with RAII semantics
pub struct GpuBuffer<T: Default + Clone + std::marker::Send = i32> {
    #[cfg(feature = "cuda")]
    ptr: Option<DevicePtr<T>>,
    #[cfg(not(feature = "cuda"))]
    ptr: Option<u64>,
    size: usize,
    #[cfg(feature = "cuda")]
    device: Arc<CudaDevice>,
    #[cfg(not(feature = "cuda"))]
    device: (),
    allocated: Arc<Mutex<u64>>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: Default + Clone + std::marker::Send> GpuBuffer<T> {
    /// Get buffer size in elements
    pub fn size(&self) -> usize {
        self.size
    }

    /// Get buffer size in bytes
    pub fn size_bytes(&self) -> usize {
        self.size * std::mem::size_of::<T>()
    }

    /// Get raw pointer
    #[cfg(feature = "cuda")]
    pub fn ptr(&self) -> Option<&DevicePtr<T>> {
        self.ptr.as_ref()
    }

    /// Get raw pointer
    #[cfg(not(feature = "cuda"))]
    pub fn ptr(&self) -> Option<u64> {
        self.ptr
    }
}

impl<T: Default + Clone + std::marker::Send> Drop for GpuBuffer<T> {
    fn drop(&mut self) {
        if let Ok(mut allocated) = self.allocated.lock() {
            *allocated = (*allocated).saturating_sub((self.size * std::mem::size_of::<T>()) as u64);
        }

        #[cfg(feature = "cuda")]
        {
            if let Some(_ptr) = self.ptr.take() {
                // Memory freed automatically when ptr is dropped
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_device_detection() {
        let _devices = GpuRuntime::detect_available_devices();
        // Just ensure it doesn't panic
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn test_gpu_runtime_creation() -> Result<()> {
        let devices = GpuRuntime::detect_available_devices()?;
        if devices.is_empty() {
            return Ok(());
        }

        let _runtime = GpuRuntime::new(devices[0])?;
        Ok(())
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn test_memory_allocation() -> Result<()> {
        let devices = GpuRuntime::detect_available_devices()?;
        if devices.is_empty() {
            return Ok(());
        }

        let runtime = GpuRuntime::new(devices[0])?;
        let initial_allocated = runtime.allocated_memory();

        let _buf: GpuBuffer<i32> = runtime.allocate(1024)?;

        let new_allocated = runtime.allocated_memory();
        assert!(new_allocated > initial_allocated);

        Ok(())
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn test_h2d_transfer() -> Result<()> {
        let devices = GpuRuntime::detect_available_devices()?;
        if devices.is_empty() {
            return Ok(());
        }

        let runtime = GpuRuntime::new(devices[0])?;
        let data = vec![1i32, 2, 3, 4, 5];

        let _gpu_buf = runtime.copy_to_device(&data)?;

        Ok(())
    }

    #[cfg(feature = "cuda")]
    #[test]
    fn test_d2h_transfer() -> Result<()> {
        let devices = GpuRuntime::detect_available_devices()?;
        if devices.is_empty() {
            return Ok(());
        }

        let runtime = GpuRuntime::new(devices[0])?;
        let data = vec![10i32, 20, 30, 40, 50];

        let gpu_buf = runtime.copy_to_device(&data)?;
        let retrieved = runtime.copy_from_device(&gpu_buf)?;

        assert_eq!(data, retrieved);

        Ok(())
    }
}
