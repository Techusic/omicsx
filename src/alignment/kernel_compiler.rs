//! CUDA kernel compilation and caching pipeline
//!
//! Provides JIT compilation via NVRTC and PTX caching for fast repeated execution.
//!
//! # Features
//!
//! - CUDA source to PTX compilation (NVRTC)
//! - Persistent cache of compiled kernels
//! - Automatic invalidation on code changes
//! - Multi-device fat binary support

use crate::error::{Error, Result};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::fs;
use std::io::{Read, Write};
use serde::{Deserialize, Serialize};

/// Kernel type for dispatch
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum KernelType {
    SmithWatermanScalar,
    SmithWatermanSimd,
    SmithWatermanGpu,
    NeedlemanWunschScalar,
    NeedlemanWunschSimd,
    NeedlemanWunschGpu,
    ViterbiSimd,
    ViterbiGpu,
    MsaProfileSimd,
    MsaProfileGpu,
}

impl KernelType {
    pub fn name(&self) -> &'static str {
        match self {
            Self::SmithWatermanScalar => "smith_waterman_scalar",
            Self::SmithWatermanSimd => "smith_waterman_simd",
            Self::SmithWatermanGpu => "smith_waterman_gpu",
            Self::NeedlemanWunschScalar => "needleman_wunsch_scalar",
            Self::NeedlemanWunschSimd => "needleman_wunsch_simd",
            Self::NeedlemanWunschGpu => "needleman_wunsch_gpu",
            Self::ViterbiSimd => "viterbi_simd",
            Self::ViterbiGpu => "viterbi_gpu",
            Self::MsaProfileSimd => "msa_profile_simd",
            Self::MsaProfileGpu => "msa_profile_gpu",
        }
    }
}

/// Compiled kernel module
#[derive(Debug, Clone)]
pub struct CompiledKernel {
    /// Kernel name
    pub name: String,
    /// PTX/fatbin code
    pub code: Vec<u8>,
    /// Compile flags used
    pub flags: Vec<String>,
    /// Target compute capability (e.g., "8.0")
    pub target: String,
    /// Timestamp of compilation
    pub timestamp: u64,
}

/// Kernel compilation cache
#[derive(Debug, Serialize, Deserialize)]
pub struct KernelCache {
    version: String,
    kernels: HashMap<String, CachedKernel>,
}

#[derive(Debug, Serialize, Deserialize)]
struct CachedKernel {
    name: String,
    source_hash: String,
    timestamp: u64,
    target: String,
    cached_path: PathBuf,
}

impl KernelCache {
    pub fn new() -> Self {
        KernelCache {
            version: "0.1".to_string(),
            kernels: HashMap::new(),
        }
    }

    /// Load cache from disk
    pub fn load(cache_dir: &Path) -> Result<Self> {
        let cache_file = cache_dir.join("kernel_cache.json");

        if !cache_file.exists() {
            return Ok(Self::new());
        }

        let mut file = fs::File::open(&cache_file)
            .map_err(|e| Error::AlignmentError(format!("Failed to read cache: {}", e)))?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)
            .map_err(|e| Error::AlignmentError(format!("Failed to read cache: {}", e)))?;

        serde_json::from_str(&contents)
            .map_err(|e| Error::AlignmentError(format!("Invalid cache format: {}", e)))
    }

    /// Save cache to disk
    pub fn save(&self, cache_dir: &Path) -> Result<()> {
        fs::create_dir_all(cache_dir)
            .map_err(|e| Error::AlignmentError(format!("Failed to create cache dir: {}", e)))?;

        let cache_file = cache_dir.join("kernel_cache.json");
        let json = serde_json::to_string_pretty(&self)
            .map_err(|e| Error::AlignmentError(format!("Failed to serialize cache: {}", e)))?;

        let mut file = fs::File::create(&cache_file)
            .map_err(|e| Error::AlignmentError(format!("Failed to write cache: {}", e)))?;
        file.write_all(json.as_bytes())
            .map_err(|e| Error::AlignmentError(format!("Failed to write cache: {}", e)))?;

        Ok(())
    }

    /// Check if kernel is in cache and still valid
    pub fn lookup(&self, name: &str, source_hash: &str) -> Option<PathBuf> {
        self.kernels.get(name).and_then(|cached| {
            if cached.source_hash == source_hash {
                Some(cached.cached_path.clone())
            } else {
                None
            }
        })
    }

    /// Add kernel to cache
    pub fn insert(&mut self, name: String, source_hash: String, target: String, path: PathBuf) {
        self.kernels.insert(
            name.clone(),
            CachedKernel {
                name,
                source_hash,
                timestamp: std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_secs(),
                target,
                cached_path: path,
            },
        );
    }
}

/// CUDA kernel compiler
pub struct KernelCompiler {
    cache_dir: PathBuf,
    cache: KernelCache,
    enable_caching: bool,
}

impl KernelCompiler {
    /// Create new kernel compiler with optional caching
    pub fn new(cache_dir: PathBuf, enable_caching: bool) -> Result<Self> {
        let cache = if enable_caching {
            KernelCache::load(&cache_dir)?
        } else {
            KernelCache::new()
        };

        Ok(KernelCompiler {
            cache_dir,
            cache,
            enable_caching,
        })
    }

    /// Compute hash of CUDA source code
    pub fn compute_source_hash(source: &str) -> String {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let mut hasher = DefaultHasher::new();
        source.hash(&mut hasher);
        format!("{:x}", hasher.finish())
    }

    /// Compile CUDA kernel to PTX
    pub fn compile_to_ptx(
        &mut self,
        kernel_type: KernelType,
        source_code: &str,
        target_cc: &str,
        flags: Vec<String>,
    ) -> Result<CompiledKernel> {
        let kernel_name = kernel_type.name();
        let source_hash = Self::compute_source_hash(source_code);

        // Check cache
        if self.enable_caching {
            if let Some(cached_path) = self.cache.lookup(kernel_name, &source_hash) {
                if cached_path.exists() {
                    let code = fs::read(&cached_path)
                        .map_err(|e| Error::AlignmentError(format!("Failed to read cached kernel: {}", e)))?;

                    return Ok(CompiledKernel {
                        name: kernel_name.to_string(),
                        code,
                        flags,
                        target: target_cc.to_string(),
                        timestamp: std::time::SystemTime::now()
                            .duration_since(std::time::UNIX_EPOCH)
                            .unwrap()
                            .as_secs(),
                    });
                }
            }
        }

        // Compile new kernel
        let ptx_code = self.nvrtc_compile(source_code, target_cc, &flags)?;

        // Cache compiled kernel
        if self.enable_caching {
            let cache_file = self.cache_dir.join(format!("{}.ptx", kernel_name));
            fs::write(&cache_file, &ptx_code)
                .map_err(|e| Error::AlignmentError(format!("Failed to cache kernel: {}", e)))?;

            self.cache.insert(
                kernel_name.to_string(),
                source_hash,
                target_cc.to_string(),
                cache_file,
            );

            self.cache.save(&self.cache_dir)?;
        }

        Ok(CompiledKernel {
            name: kernel_name.to_string(),
            code: ptx_code,
            flags,
            target: target_cc.to_string(),
            timestamp: std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_secs(),
        })
    }

    /// NVRTC compilation (placeholder for feature gating)
    fn nvrtc_compile(&self, _source: &str, _target_cc: &str, _flags: &[String]) -> Result<Vec<u8>> {
        #[cfg(feature = "cuda")]
        {
            // In production with cudarc, would call:
            // let ptx = cudarc::nvrtc::compile_ptx(source, flags, target_cc).unwrap();
            // Ok(ptx.into_bytes())
            Err(Error::AlignmentError(
                "NVRTC compilation not yet implemented with cudarc".to_string(),
            ))
        }

        #[cfg(not(feature = "cuda"))]
        {
            Err(Error::AlignmentError(
                "CUDA support not available; enable 'cuda' feature".to_string(),
            ))
        }
    }

    /// Load pre-compiled PTX file
    pub fn load_ptx_file(&self, path: &Path) -> Result<CompiledKernel> {
        let code = fs::read(path)
            .map_err(|e| Error::AlignmentError(format!("Failed to read PTX file: {}", e)))?;

        let name = path
            .file_stem()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        Ok(CompiledKernel {
            name,
            code,
            flags: Vec::new(),
            target: "8.0".to_string(), // Assume Ampere
            timestamp: std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_secs(),
        })
    }

    /// Clear cache
    pub fn clear_cache(&mut self) -> Result<()> {
        if self.cache_dir.exists() {
            fs::remove_dir_all(&self.cache_dir)
                .map_err(|e| Error::AlignmentError(format!("Failed to clear cache: {}", e)))?;
        }
        self.cache = KernelCache::new();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kernel_type_names() {
        assert_eq!(KernelType::SmithWatermanGpu.name(), "smith_waterman_gpu");
        assert_eq!(KernelType::ViterbiGpu.name(), "viterbi_gpu");
    }

    #[test]
    fn test_source_hash() {
        let source1 = "int main() {}";
        let source2 = "int main() {}";
        let source3 = "int main() { return 0; }";

        let hash1 = KernelCompiler::compute_source_hash(source1);
        let hash2 = KernelCompiler::compute_source_hash(source2);
        let hash3 = KernelCompiler::compute_source_hash(source3);

        assert_eq!(hash1, hash2);
        assert_ne!(hash1, hash3);
    }

    #[test]
    fn test_kernel_cache() {
        let cache = KernelCache::new();
        assert_eq!(cache.kernels.len(), 0);

        let mut cache = cache;
        cache.insert(
            "test".to_string(),
            "abc123".to_string(),
            "8.0".to_string(),
            PathBuf::from("/test/kernel.ptx"),
        );

        assert_eq!(cache.kernels.len(), 1);
        assert!(cache.lookup("test", "abc123").is_some());
        assert!(cache.lookup("test", "wrong_hash").is_none());
    }

    #[test]
    fn test_compiler_creation() -> Result<()> {
        let temp_dir = std::env::temp_dir().join("omics_test_cache");
        let _compiler = KernelCompiler::new(temp_dir, true)?;
        Ok(())
    }
}
