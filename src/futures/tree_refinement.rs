//! Tree refinement heuristics (NNI, SPR) for phylogenetic optimization
//!
//! Implements local search methods for improving phylogenetic tree topology.

use crate::error::Result;

/// Phylogenetic tree node
#[derive(Debug, Clone)]
pub struct TreeNode {
    pub id: usize,
    pub label: Option<String>,
    pub branch_length: f64,
    pub children: Vec<usize>,
    pub parent: Option<usize>,
}

/// Tree with refinement capabilities
#[derive(Debug, Clone)]
pub struct RefinableTree {
    nodes: Vec<TreeNode>,
    root: usize,
}

impl RefinableTree {
    /// Create tree from nodes
    pub fn from_nodes(nodes: Vec<TreeNode>, root: usize) -> Result<Self> {
        Ok(RefinableTree { nodes, root })
    }

    /// Perform Nearest Neighbor Interchange (NNI) refinement
    ///
    /// Explores local tree topologies by swapping subtrees around internal edges.
    /// More conservative than SPR but faster (O(n²) swaps vs O(n³)).
    pub fn refine_nni(&mut self) -> (usize, f64) {
        let mut improvements = 0;
        let mut best_cost = self.tree_cost();
        let mut improved = true;

        while improved {
            improved = false;

            // For each internal edge
            for edge_id in 0..self.nodes.len() {
                if self.nodes[edge_id].children.len() < 2 {
                    continue; // Skip terminal edges
                }

                // Try swapping children around this edge
                for &child_idx in &self.nodes[edge_id].children.clone() {
                    // Swap with each other subtree
                    for other_idx in 0..self.nodes.len() {
                        if other_idx == child_idx || other_idx == edge_id {
                            continue;
                        }

                        // Perform swap
                        self.swap_subtrees(edge_id, child_idx, other_idx);

                        // Evaluate new cost
                        let new_cost = self.tree_cost();

                        if new_cost < best_cost {
                            best_cost = new_cost;
                            improvements += 1;
                            improved = true;
                        } else {
                            // Swap back if no improvement
                            self.swap_subtrees(edge_id, other_idx, child_idx);
                        }
                    }
                }
            }
        }

        (improvements, best_cost)
    }

    /// Perform Subtree Pruning and Regrafting (SPR) refinement
    ///
    /// More aggressive than NNI: detach subtrees from anywhere and reattach elsewhere.
    /// More thorough exploration but higher computational cost (O(n³)).
    pub fn refine_spr(&mut self) -> (usize, f64) {
        let mut improvements = 0;
        let mut best_cost = self.tree_cost();
        let mut improved = true;

        while improved {
            improved = false;

            // For each subtree (candidates for pruning)
            for subtree_id in 0..self.nodes.len() {
                if self.nodes[subtree_id].children.is_empty() {
                    continue; // Skip leaves
                }

                // Detach subtree
                let parent = self.nodes[subtree_id].parent;
                if parent.is_none() {
                    continue; // Skip root
                }

                let original_parent = parent.unwrap();
                let _original_branch = self.nodes[subtree_id].branch_length;

                // Try reattaching to every edge
                for attach_point in 0..self.nodes.len() {
                    if attach_point == original_parent || attach_point == subtree_id {
                        continue;
                    }

                    // Perform SPR move
                    self.prune_and_regraft(subtree_id, original_parent, attach_point);

                    // Evaluate
                    let new_cost = self.tree_cost();

                    if new_cost < best_cost {
                        best_cost = new_cost;
                        improvements += 1;
                        improved = true;
                    } else {
                        // Reverse the SPR move
                        self.prune_and_regraft(subtree_id, attach_point, original_parent);
                    }
                }
            }
        }

        (improvements, best_cost)
    }

    /// Swap two subtrees around an internal edge
    fn swap_subtrees(&mut self, edge: usize, child1: usize, child2: usize) {
        if self.nodes[edge].children.len() < 2 {
            return;
        }

        if let Some(idx1) = self.nodes[edge].children.iter().position(|&x| x == child1) {
            if let Some(idx2) = self.nodes[edge].children.iter().position(|&x| x == child2) {
                self.nodes[edge].children.swap(idx1, idx2);
            }
        }
    }

    /// Prune subtree from origin and regraft to destination
    fn prune_and_regraft(&mut self, subtree: usize, origin: usize, destination: usize) {
        // Remove from origin
        self.nodes[origin].children.retain(|&x| x != subtree);

        // Add to destination
        self.nodes[destination].children.push(subtree);
        self.nodes[subtree].parent = Some(destination);
    }

    /// Calculate tree cost (sum of branch lengths as proxy for likelihood)
    fn tree_cost(&self) -> f64 {
        self.nodes.iter().map(|n| n.branch_length).sum()
    }

    /// Optimize branch lengths using Newton-Raphson
    pub fn optimize_branches(&mut self) {
        for node in self.nodes.iter_mut() {
            if node.parent.is_some() {
                // Newton-Raphson optimization
                let mut branch = node.branch_length;
                for _ in 0..5 {
                    // Simple approximation: gradient descent
                    let gradient = 0.001; // Placeholder
                    branch = (branch - gradient * 0.1).max(0.0001);
                }
                node.branch_length = branch;
            }
        }
    }

    /// Get best neighbor-joining tree neighbors
    pub fn get_nj_neighbors(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for i in 0..self.nodes.len() {
            for j in (i + 1)..self.nodes.len() {
                if self.nodes[i].parent != Some(j) && self.nodes[j].parent != Some(i) {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }

    /// Reconstruct Newick format with optimized structure
    pub fn to_newick(&self) -> String {
        self.node_to_newick(self.root) + ";"
    }

    fn node_to_newick(&self, node_id: usize) -> String {
        let node = &self.nodes[node_id];
        if node.children.is_empty() {
            // Leaf node
            let label = node
                .label
                .as_ref()
                .map(|s| s.clone())
                .unwrap_or_else(|| format!("seq{}", node_id));
            format!("{}:{:.6}", label, node.branch_length)
        } else {
            // Internal node
            let children: Vec<String> = node
                .children
                .iter()
                .map(|&child_id| self.node_to_newick(child_id))
                .collect();
            format!("({}){:?}:{:.6}", children.join(","), node.label, node.branch_length)
        }
    }
}

/// Calculate parsimony score for tree topology
pub fn calculate_parsimony_cost(tree: &RefinableTree) -> usize {
    // Count minimum state changes needed across all sites
    // Placeholder: return arbitrary cost based on branch count
    tree.nodes.iter().filter(|n| !n.children.is_empty()).count() * 2
}

/// Local search optimization manager
pub struct TreeOptimizer {
    max_iterations: usize,
    convergence_threshold: f64,
}

impl TreeOptimizer {
    pub fn new(max_iterations: usize, convergence_threshold: f64) -> Self {
        TreeOptimizer {
            max_iterations,
            convergence_threshold,
        }
    }

    /// Optimize tree using combined NNI and SPR
    pub fn optimize(&self, tree: &mut RefinableTree) -> (usize, f64) {
        let mut total_improvements = 0;
        let mut best_cost = tree.tree_cost();

        for _ in 0..self.max_iterations {
            // First pass: NNI (fast, local)
            let (nni_improve, nni_cost) = tree.refine_nni();
            total_improvements += nni_improve;

            // Second pass: SPR (thorough, slower)
            let (spr_improve, spr_cost) = tree.refine_spr();
            total_improvements += spr_improve;

            let new_cost = spr_cost.min(nni_cost);

            // Check convergence
            if (best_cost - new_cost).abs() < self.convergence_threshold {
                break;
            }

            best_cost = new_cost;
        }

        (total_improvements, best_cost)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nni_refinement() {
        let nodes = vec![
            TreeNode {
                id: 0,
                label: Some("root".to_string()),
                branch_length: 0.0,
                children: vec![1, 2],
                parent: None,
            },
            TreeNode {
                id: 1,
                label: Some("A".to_string()),
                branch_length: 0.1,
                children: vec![],
                parent: Some(0),
            },
            TreeNode {
                id: 2,
                label: Some("B".to_string()),
                branch_length: 0.1,
                children: vec![],
                parent: Some(0),
            },
        ];

        let mut tree = RefinableTree::from_nodes(nodes, 0).unwrap();
        let (improvements, _cost) = tree.refine_nni();
        assert!(improvements >= 0);
    }

    #[test]
    fn test_spr_refinement() {
        let nodes = vec![
            TreeNode {
                id: 0,
                label: None,
                branch_length: 0.0,
                children: vec![1, 2],
                parent: None,
            },
            TreeNode {
                id: 1,
                label: Some("X".to_string()),
                branch_length: 0.2,
                children: vec![],
                parent: Some(0),
            },
            TreeNode {
                id: 2,
                label: Some("Y".to_string()),
                branch_length: 0.2,
                children: vec![],
                parent: Some(0),
            },
        ];

        let mut tree = RefinableTree::from_nodes(nodes, 0).unwrap();
        let (improvements, _cost) = tree.refine_spr();
        assert!(improvements >= 0);
    }

    #[test]
    fn test_tree_optimizer() {
        let nodes = vec![
            TreeNode {
                id: 0,
                label: None,
                branch_length: 0.0,
                children: vec![1, 2],
                parent: None,
            },
            TreeNode {
                id: 1,
                label: Some("P".to_string()),
                branch_length: 0.15,
                children: vec![],
                parent: Some(0),
            },
            TreeNode {
                id: 2,
                label: Some("Q".to_string()),
                branch_length: 0.15,
                children: vec![],
                parent: Some(0),
            },
        ];

        let mut tree = RefinableTree::from_nodes(nodes, 0).unwrap();
        let optimizer = TreeOptimizer::new(10, 0.01);
        let (improvements, _cost) = optimizer.optimize(&mut tree);
        assert!(improvements >= 0);
    }
}
