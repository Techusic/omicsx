//! GPU Performance Benchmarks
//!
//! Compares performance across scalar, SIMD, and GPU implementations.
//! Run with: cargo bench --bench gpu_benchmarks -- --verbose

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use omics_simd::alignment::{GpuDispatcher, AlignmentStrategy};
use omics_simd::alignment::gpu_dispatcher::GpuDispatcherStrategy;

fn gpu_strategy_selection_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("gpu_strategy_selection");
    
    for size in [100, 1000, 10000, 100000].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| {
                GpuDispatcherStrategy::select_strategy(black_box(size), black_box(size), true, None)
            });
        });
    }
    group.finish();
}

fn gpu_memory_estimation_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("gpu_memory_estimation");
    
    for size in [100, 1000, 10000, 100000].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| {
                GpuDispatcherStrategy::estimate_gpu_memory(black_box(size), black_box(size))
            });
        });
    }
    group.finish();
}

fn gpu_dispatcher_initialization(c: &mut Criterion) {
    c.bench_function("gpu_dispatcher_new", |b| {
        b.iter(|| {
            GpuDispatcher::new()
        });
    });
}

fn gpu_fitness_check_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("gpu_fitness_check");
    let available_memory = 8 * 1024 * 1024 * 1024; // 8GB
    
    for size in [100, 1000, 10000, 100000].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| {
                GpuDispatcherStrategy::fits_in_gpu_memory(
                    black_box(size),
                    black_box(size),
                    black_box(available_memory),
                )
            });
        });
    }
    group.finish();
}

fn cuda_kernel_timing_simulation(c: &mut Criterion) {
    // Simulated CUDA kernel timings (in production, would measure actual kernel execution)
    c.bench_function("cuda_kernel_overhead", |b| {
        b.iter(|| {
            // Simulate kernel launch overhead (~microseconds)
            #[cfg(feature = "cuda")]
            {
                use omics_simd::alignment::kernel::cuda::CudaAlignmentKernel;
                let _kernel = CudaAlignmentKernel::new();
            }
        });
    });
}

fn hip_kernel_timing_simulation(c: &mut Criterion) {
    c.bench_function("hip_kernel_overhead", |b| {
        b.iter(|| {
            #[cfg(feature = "hip")]
            {
                use omics_simd::alignment::kernel::hip::HipAlignmentKernel;
                let _kernel = HipAlignmentKernel::new();
            }
        });
    });
}

fn vulkan_kernel_timing_simulation(c: &mut Criterion) {
    c.bench_function("vulkan_kernel_overhead", |b| {
        b.iter(|| {
            #[cfg(feature = "vulkan")]
            {
                use omics_simd::alignment::kernel::vulkan::VulkanComputeKernel;
                let _kernel = VulkanComputeKernel::new();
            }
        });
    });
}

fn speedup_estimation_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("speedup_estimation");
    
    let strategies = vec![
        ("Scalar", AlignmentStrategy::Scalar),
        ("SIMD", AlignmentStrategy::Simd),
        ("Banded", AlignmentStrategy::Banded),
        ("GPU Full", AlignmentStrategy::GpuFull),
        ("GPU Tiled", AlignmentStrategy::GpuTiled),
    ];
    
    for (name, strategy) in strategies {
        group.bench_with_input(BenchmarkId::from_parameter(name), &strategy, |b, &strategy| {
            b.iter(|| {
                GpuDispatcherStrategy::gpu_speedup_factor(black_box(strategy))
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    gpu_dispatcher_initialization,
    gpu_strategy_selection_benchmark,
    gpu_memory_estimation_benchmark,
    gpu_fitness_check_benchmark,
    cuda_kernel_timing_simulation,
    hip_kernel_timing_simulation,
    vulkan_kernel_timing_simulation,
    speedup_estimation_benchmark,
);
criterion_main!(benches);
