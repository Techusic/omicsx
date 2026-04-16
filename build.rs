// Build script for OMICSX - CUDA detection and configuration
// 
// CUDA Support Information:
// ========================
// By default, OMICSX builds without GPU support to ensure maximum compatibility.
// 
// To enable CUDA GPU acceleration:
//   1. Ensure you have CUDA 12.5 or earlier installed
//   2. Build with: cargo build --release --features cuda
//
// Note: cudarc (the CUDA Rust library) doesn't support CUDA 13.2+ yet.
// If you have CUDA 13.x:
//   - Option 1: Downgrade to CUDA 12.5.x for full GPU support
//   - Option 2: Use scalar/NEON backends (default, still fast with SIMD)

use std::path::Path;
use std::process::Command;

fn main() {
    // Check if CUDA is available (informational only)
    if let Ok(cuda_path) = std::env::var("CUDA_PATH") {
        if Path::new(&cuda_path).exists() {
            if let Some(version) = get_cuda_version(&cuda_path) {
                eprintln!("[OMICSX build] CUDA {} detected", version);

                // Check for known incompatibilities
                if version.starts_with("13") && !version.starts_with("13.0") && !version.starts_with("13.1") {
                    eprintln!("[OMICSX build] ⚠️  CUDA {} is not yet supported by cudarc", version);
                    eprintln!("[OMICSX build]     GPU features will be disabled automatically");
                    eprintln!("[OMICSX build]     To enable GPU support, downgrade to CUDA 12.5.x or use NEON backend");
                }
            }
        }
    }

    // Rerun if CUDA environment changes
    println!("cargo:rerun-if-env-changed=CUDA_PATH");
    println!("cargo:rerun-if-env-changed=CUDA_ROOT");
}

/// Get CUDA version from nvcc compiler
fn get_cuda_version(cuda_path: &str) -> Option<String> {
    let nvcc = if cfg!(windows) {
        Path::new(cuda_path).join("bin").join("nvcc.exe")
    } else {
        Path::new(cuda_path).join("bin").join("nvcc")
    };

    if !nvcc.exists() {
        return None;
    }

    let output = Command::new(&nvcc)
        .arg("--version")
        .output()
        .ok()?;

    let stdout = String::from_utf8(output.stdout).ok()?;
    
    for line in stdout.lines() {
        if line.contains("release") {
            if let Some(version_str) = line.split("release ").nth(1) {
                if let Some(version) = version_str.trim().split_whitespace().next() {
                    return Some(version.to_string());
                }
            }
        }
    }

    None
}
