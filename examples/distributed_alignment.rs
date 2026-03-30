/// Distributed Multi-Node Alignment Demonstration
///
/// This example demonstrates how to use the distributed coordinator
/// to parallelize sequence alignment across multiple nodes with automatic
/// work distribution and load balancing.

use omicsx::futures::DistributedCoordinator;
use omicsx::futures::AlignmentTask;
use omicsx::protein::Protein;
use omicsx::scoring::ScoringMatrix;
use std::thread;
use std::time::Duration;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== OMICSX Distributed Multi-Node Alignment Demo ===\n");

    // Create a distributed coordinator
    let coordinator = DistributedCoordinator::new();
    println!("✓ Created distributed coordinator\n");

    // Register multiple nodes in the cluster
    let node_ids: Vec<_> = (0..4)
        .map(|i| {
            let node_id = coordinator.register_node().unwrap();
            println!("  Registered Node {}: {:?}", i, node_id);
            node_id
        })
        .collect();

    println!("\n{} nodes registered\n", node_ids.len());

    // Create a batch of alignment tasks
    let seq1 = "MTEITAAMAIDQEAIDAAAGHDKLSLVGANKEDIAAANDLGAWILGHPDNHEAVGVQVKVKALPDAQFEVVHSLAKWKRQTLG";
    let seq2 = "MGQVGRQLGHHDLSAGNEEPGAADLSAGCQVGRQLGHHDLSAGNEEPGAADLSAGAHGDVPPEDQVGQQVGHDLSVGGDHLLEHSQGAQ";

    let tasks: Vec<AlignmentTask> = (0..20)
        .map(|i| {
            // Vary sequences for realistic test
            let query_seq = if i % 2 == 0 { seq1 } else { seq2 };
            let subject_seq = if i % 3 == 0 { seq2 } else { seq1 };

            AlignmentTask {
                task_id: i,
                query: Protein::from_string(query_seq)
                    .unwrap()
                    .with_id(format!("query_{}", i)),
                subject: Protein::from_string(subject_seq)
                    .unwrap()
                    .with_id(format!("subject_{}", i)),
                matrix: ScoringMatrix::default(),
            }
        })
        .collect();

    // Submit the batch
    let submitted = coordinator.submit_batch(tasks)?;
    println!("Submitted {} alignment tasks\n", submitted);

    // Show initial stats
    let initial_stats = coordinator.get_stats()?;
    println!("Initial Distribution Stats:");
    println!("  Total nodes: {}", initial_stats.total_nodes);
    println!("  Pending tasks: {}\n", initial_stats.pending_tasks);

    // Simulate node workers pulling tasks (work-stealing distribution)
    println!("Simulating distributed work processing...\n");

    for (idx, node_id) in node_ids.iter().enumerate() {
        println!("Node {}: Starting work...", idx);

        // Simulate this node stealing work
        loop {
            match coordinator.get_task(*node_id) {
                Ok(Some(task)) => {
                    // Simulate alignment processing
                    thread::sleep(Duration::from_millis(10));

                    // Record result (simulated)
                    let result = omicsx::futures::AlignmentResultRecord {
                        task_id: task.task_id,
                        node_id: *node_id,
                        score: 100 + (task.task_id as i32 * 5),
                        identity: 0.85 + (task.task_id as f32 * 0.01),
                        gaps: task.task_id,
                        query_coverage: 0.95,
                    };

                    coordinator.record_result(result)?;
                    println!("  Node {}: Completed task {}", idx, task.task_id);
                }
                Ok(None) => {
                    println!("  Node {}: No more tasks\n", idx);
                    break;
                }
                Err(e) => {
                    eprintln!("  Node {}: Error: {}\n", idx, e);
                    break;
                }
            }
        }
    }

    // Show final statistics
    println!("Distribution Complete!\n");

    let final_stats = coordinator.get_stats()?;
    println!("Final Distribution Stats:");
    println!("  Total nodes: {}", final_stats.total_nodes);
    println!("  Online nodes: {}", final_stats.online_nodes);
    println!("  Completed tasks: {}", final_stats.completed_tasks);
    println!("  Pending tasks: {}", final_stats.pending_tasks);
    println!("  Total compute time: {} ms", final_stats.total_compute_time_ms);
    println!("  Average time per task: {} ms\n", final_stats.average_time_per_task_ms);

    // Show per-node statistics
    println!("Node Statistics:");
    let node_stats = coordinator.get_node_stats()?;
    for stat in node_stats {
        println!("  Node {:?}:", stat.node_id);
        println!("    Tasks: {}", stat.task_count);
        println!("    Completed: {}", stat.completed_tasks);
        println!("    Status: {:?}", stat.status);
    }

    println!("\n// Key Features Demonstrated:");
    println!("// 1. Multi-node registration and management");
    println!("// 2. Batch task submission");
    println!("// 3. Work-stealing load balancing");
    println!("// 4. Result aggregation");
    println!("// 5. Statistical tracking and reporting");

    println!("\n=== Demo Complete ===");

    Ok(())
}
