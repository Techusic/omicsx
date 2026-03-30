/// Distributed multi-node alignment coordination
///
/// Provides infrastructure for parallelizing sequence alignment across multiple nodes
/// with automatic work distribution, load balancing, and result aggregation.

use crate::error::{Error, Result};
use crate::protein::Protein;
use crate::scoring::ScoringMatrix;
use std::collections::{HashMap, VecDeque};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Unique identifier for nodes in the distributed system
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct NodeId(pub usize);

/// Alignment task for distributed processing
#[derive(Debug, Clone)]
pub struct AlignmentTask {
    pub task_id: usize,
    pub query: Protein,
    pub subject: Protein,
    pub matrix: ScoringMatrix,
}

/// Result of a completed alignment task
#[derive(Debug, Clone)]
pub struct AlignmentResultRecord {
    pub task_id: usize,
    pub node_id: NodeId,
    pub score: i32,
    pub identity: f32,
    pub gaps: usize,
    pub query_coverage: f32,
}

/// Node statistics for load balancing
#[derive(Debug, Clone)]
pub struct NodeStats {
    pub node_id: NodeId,
    pub task_count: usize,
    pub completed_tasks: usize,
    pub total_time_ms: u128,
    pub status: NodeStatus,
}

/// Current status of a node
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeStatus {
    /// Node is ready to accept work
    Ready,
    /// Node is currently processing tasks
    Processing,
    /// Node is unreachable
    Unavailable,
    /// Node has shutdown
    Offline,
}

/// Work-stealing task queue for load balancing
pub struct TaskQueue {
    tasks: Arc<Mutex<VecDeque<AlignmentTask>>>,
    task_counter: Arc<AtomicUsize>,
}

impl TaskQueue {
    /// Create a new task queue
    pub fn new() -> Self {
        TaskQueue {
            tasks: Arc::new(Mutex::new(VecDeque::new())),
            task_counter: Arc::new(AtomicUsize::new(0)),
        }
    }

    /// Enqueue a batch of tasks
    pub fn enqueue_batch(&self, tasks: Vec<AlignmentTask>) -> Result<usize> {
        let mut queue = self.tasks.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to acquire task queue lock: {}", e))
        })?;

        let count = tasks.len();
        for task in tasks {
            queue.push_back(task);
        }

        Ok(count)
    }

    /// Dequeue a single task (work-stealing)
    pub fn dequeue(&self) -> Result<Option<AlignmentTask>> {
        let mut queue = self.tasks.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to acquire task queue lock: {}", e))
        })?;

        Ok(queue.pop_front())
    }

    /// Steal multiple tasks for batch processing
    pub fn steal_batch(&self, batch_size: usize) -> Result<Vec<AlignmentTask>> {
        let mut queue = self.tasks.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to acquire task queue lock: {}", e))
        })?;

        let mut batch = Vec::new();
        for _ in 0..batch_size {
            if let Some(task) = queue.pop_front() {
                batch.push(task);
            } else {
                break;
            }
        }

        Ok(batch)
    }

    /// Get queue length
    pub fn len(&self) -> Result<usize> {
        let queue = self.tasks.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to acquire task queue lock: {}", e))
        })?;

        Ok(queue.len())
    }

    /// Check if queue is empty
    pub fn is_empty(&self) -> Result<bool> {
        Ok(self.len()? == 0)
    }

    /// Get next task ID
    fn next_task_id(&self) -> usize {
        self.task_counter.fetch_add(1, Ordering::SeqCst)
    }
}

/// Distributed alignment coordinator
pub struct DistributedCoordinator {
    nodes: Arc<Mutex<HashMap<NodeId, NodeStats>>>,
    task_queue: TaskQueue,
    results: Arc<Mutex<Vec<AlignmentResultRecord>>>,
    next_node_id: Arc<AtomicUsize>,
}

impl DistributedCoordinator {
    /// Create a new distributed coordinator
    pub fn new() -> Self {
        DistributedCoordinator {
            nodes: Arc::new(Mutex::new(HashMap::new()),),
            task_queue: TaskQueue::new(),
            results: Arc::new(Mutex::new(Vec::new())),
            next_node_id: Arc::new(AtomicUsize::new(0)),
        }
    }

    /// Register a new node
    pub fn register_node(&self) -> Result<NodeId> {
        let node_id = NodeId(self.next_node_id.fetch_add(1, Ordering::SeqCst));

        let mut nodes = self.nodes.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to register node: {}", e))
        })?;

        nodes.insert(
            node_id,
            NodeStats {
                node_id,
                task_count: 0,
                completed_tasks: 0,
                total_time_ms: 0,
                status: NodeStatus::Ready,
            },
        );

        Ok(node_id)
    }

    /// Submit a batch of alignment tasks for distributed processing
    pub fn submit_batch(&self, tasks: Vec<AlignmentTask>) -> Result<usize> {
        self.task_queue.enqueue_batch(tasks)
    }

    /// Get next available task (work-stealing for load balancing)
    pub fn get_task(&self, node_id: NodeId) -> Result<Option<AlignmentTask>> {
        // Update node status
        self.update_node_status(node_id, NodeStatus::Processing)?;

        // Attempt to get a task
        self.task_queue.dequeue()
    }

    /// Record a completed alignment result
    pub fn record_result(&self, result: AlignmentResultRecord) -> Result<()> {
        let mut results = self.results.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to record result: {}", e))
        })?;

        results.push(result);
        Ok(())
    }

    /// Get all completed results
    pub fn get_results(&self) -> Result<Vec<AlignmentResultRecord>> {
        let results = self.results.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to retrieve results: {}", e))
        })?;

        Ok(results.clone())
    }

    /// Get statistics for all nodes
    pub fn get_node_stats(&self) -> Result<Vec<NodeStats>> {
        let nodes = self.nodes.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to retrieve node stats: {}", e))
        })?;

        Ok(nodes.values().cloned().collect())
    }

    /// Update node status
    fn update_node_status(&self, node_id: NodeId, status: NodeStatus) -> Result<()> {
        let mut nodes = self.nodes.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to update node status: {}", e))
        })?;

        if let Some(stats) = nodes.get_mut(&node_id) {
            stats.status = status;
        }

        Ok(())
    }

    /// Get remaining task count
    pub fn pending_tasks(&self) -> Result<usize> {
        self.task_queue.len()
    }

    /// Check if all work is complete
    pub fn is_complete(&self) -> Result<bool> {
        self.task_queue.is_empty()
    }

    /// Get coordination statistics
    pub fn get_stats(&self) -> Result<DistributionStats> {
        let nodes = self.nodes.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to get stats: {}", e))
        })?;

        let results = self.results.lock().map_err(|e| {
            Error::AlignmentError(format!("Failed to get results count: {}", e))
        })?;

        let pending = self.pending_tasks()?;
        let completed = results.len();
        let online_nodes = nodes.values().filter(|s| s.status != NodeStatus::Offline).count();

        let total_time_ms = nodes.values().map(|s| s.total_time_ms).sum();

        Ok(DistributionStats {
            total_nodes: nodes.len(),
            online_nodes,
            pending_tasks: pending,
            completed_tasks: completed,
            total_compute_time_ms: total_time_ms,
            average_time_per_task_ms: if completed > 0 {
                total_time_ms / completed as u128
            } else {
                0
            },
        })
    }
}

impl Default for DistributedCoordinator {
    fn default() -> Self {
        Self::new()
    }
}

/// Statistics about the distributed alignment run
#[derive(Debug, Clone)]
pub struct DistributionStats {
    pub total_nodes: usize,
    pub online_nodes: usize,
    pub pending_tasks: usize,
    pub completed_tasks: usize,
    pub total_compute_time_ms: u128,
    pub average_time_per_task_ms: u128,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::Protein;
    use crate::scoring::ScoringMatrix;

    #[test]
    fn test_coordinator_creation() {
        let coordinator = DistributedCoordinator::new();
        assert!(coordinator.pending_tasks().is_ok());
    }

    #[test]
    fn test_node_registration() -> Result<()> {
        let coordinator = DistributedCoordinator::new();
        let node_id = coordinator.register_node()?;
        assert_eq!(node_id, NodeId(0));

        let node_id2 = coordinator.register_node()?;
        assert_eq!(node_id2, NodeId(1));

        let stats = coordinator.get_node_stats()?;
        assert_eq!(stats.len(), 2);

        Ok(())
    }

    #[test]
    fn test_task_queue() -> Result<()> {
        let queue = TaskQueue::new();
        assert!(queue.is_empty()?);

        let task = AlignmentTask {
            task_id: 0,
            query: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap().with_id("test1".to_string()),
            subject: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap().with_id("test2".to_string()),
            matrix: ScoringMatrix::default(),
        };

        queue.enqueue_batch(vec![task.clone()])?;
        assert!(!queue.is_empty()?);

        let retrieved = queue.dequeue()?;
        assert!(retrieved.is_some());
        assert!(queue.is_empty()?);

        Ok(())
    }

    #[test]
    fn test_batch_submission() -> Result<()> {
        let coordinator = DistributedCoordinator::new();
        let tasks: Vec<AlignmentTask> = (0..5)
            .map(|i| AlignmentTask {
                task_id: i,
                query: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap()
                    .with_id(format!("query_{}", i)),
                subject: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap()
                    .with_id(format!("subject_{}", i)),
                matrix: ScoringMatrix::default(),
            })
            .collect();

        let count = coordinator.submit_batch(tasks)?;
        assert_eq!(count, 5);
        assert_eq!(coordinator.pending_tasks()?, 5);

        Ok(())
    }

    #[test]
    fn test_work_stealing() -> Result<()> {
        let coordinator = DistributedCoordinator::new();
        let node_id = coordinator.register_node()?;

        let tasks: Vec<AlignmentTask> = (0..3)
            .map(|i| AlignmentTask {
                task_id: i,
                query: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap()
                    .with_id(format!("query_{}", i)),
                subject: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap()
                    .with_id(format!("subject_{}", i)),
                matrix: ScoringMatrix::default(),
            })
            .collect();

        coordinator.submit_batch(tasks)?;

        let task1 = coordinator.get_task(node_id)?;
        assert!(task1.is_some());

        let task2 = coordinator.get_task(node_id)?;
        assert!(task2.is_some());

        let task3 = coordinator.get_task(node_id)?;
        assert!(task3.is_some());

        assert!(coordinator.is_complete()?);

        Ok(())
    }

    #[test]
    fn test_result_recording() -> Result<()> {
        let coordinator = DistributedCoordinator::new();

        let result = AlignmentResultRecord {
            task_id: 0,
            node_id: NodeId(0),
            score: 100,
            identity: 0.95,
            gaps: 2,
            query_coverage: 0.98,
        };

        coordinator.record_result(result)?;
        let results = coordinator.get_results()?;
        assert_eq!(results.len(), 1);

        Ok(())
    }

    #[test]
    fn test_distribution_stats() -> Result<()> {
        let coordinator = DistributedCoordinator::new();
        coordinator.register_node()?;
        coordinator.register_node()?;

        let tasks: Vec<AlignmentTask> = (0..10)
            .map(|i| AlignmentTask {
                task_id: i,
                query: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap()
                    .with_id(format!("query_{}", i)),
                subject: Protein::from_string("ACDEFGHIKLMNPQRSTVWY").unwrap()
                    .with_id(format!("subject_{}", i)),
                matrix: ScoringMatrix::default(),
            })
            .collect();

        coordinator.submit_batch(tasks)?;

        let stats = coordinator.get_stats()?;
        assert_eq!(stats.total_nodes, 2);
        assert_eq!(stats.pending_tasks, 10);
        assert_eq!(stats.completed_tasks, 0);

        Ok(())
    }

    #[test]
    fn test_node_status_tracking() -> Result<()> {
        let coordinator = DistributedCoordinator::new();
        let node_id = coordinator.register_node()?;

        let stats = coordinator.get_node_stats()?;
        assert_eq!(stats[0].status, NodeStatus::Ready);

        Ok(())
    }
}
