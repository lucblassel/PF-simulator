use std::fmt::Write;

use crate::rng::Rng;
pub use crate::rng::sample_seed;

#[derive(Clone, Debug)]
pub struct Edge {
    pub a: usize,
    pub b: usize,
    pub length: f64,
}

/// Unrooted binary tree. Tip IDs 0..tips are labeled T1..Tn in Newick.
#[derive(Debug)]
pub struct Tree {
    tips: usize,
    edges: Vec<Edge>,
}

impl Tree {
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Uniform over unrooted labeled binary topologies, with independent
    /// exponential branch lengths measured in substitutions per site.
    pub fn uniform(tips: usize, branch_mean: f64, seed: u64) -> Result<Self, String> {
        if !(3..=10_000).contains(&tips)
            || !branch_mean.is_finite()
            || !(1e-8..=100.0).contains(&branch_mean)
        {
            return Err("Need 3 <= tips <= 10000 and a finite branch mean in [1e-8, 100]".into());
        }

        let mut rng = Rng::new(seed);
        let mut edges = Vec::with_capacity(2 * tips - 3);
        for tip in 0..3 {
            edges.push(Edge {
                a: tip,
                b: tips,
                length: 0.0,
            });
        }

        // Deleting the last tip uniquely recovers its predecessor and insertion
        // edge. Choosing edges uniformly therefore yields a uniform topology.
        for tip in 3..tips {
            let index = rng.index(edges.len());
            let old = edges[index].clone();
            let internal = tips + tip - 2;
            edges[index] = Edge {
                a: old.a,
                b: internal,
                length: 0.0,
            };
            edges.push(Edge {
                a: internal,
                b: old.b,
                length: 0.0,
            });
            edges.push(Edge {
                a: internal,
                b: tip,
                length: 0.0,
            });
        }

        for edge in &mut edges {
            edge.length = rng.exponential(branch_mean);
        }
        Ok(Self { tips, edges })
    }

    pub fn newick(&self) -> String {
        let mut neighbors = vec![Vec::new(); 2 * self.tips - 2];
        for edge in &self.edges {
            neighbors[edge.a].push((edge.b, edge.length));
            neighbors[edge.b].push((edge.a, edge.length));
        }

        // Iterative traversal avoids overflowing Rust's stack for large trees.
        enum Task {
            Node(usize, usize),
            Length(f64),
            Text(&'static str),
        }
        let mut tasks = vec![Task::Text(";\n"), Task::Node(self.tips, usize::MAX)];
        let mut result = String::new();
        while let Some(task) = tasks.pop() {
            match task {
                Task::Text(text) => result.push_str(text),
                Task::Length(length) => write!(&mut result, ":{length:.17e}").unwrap(),
                Task::Node(node, _) if node < self.tips => {
                    write!(&mut result, "T{}", node + 1).unwrap();
                }
                Task::Node(node, parent) => {
                    result.push('(');
                    tasks.push(Task::Text(")"));
                    let children: Vec<_> = neighbors[node]
                        .iter()
                        .filter(|(child, _)| *child != parent)
                        .collect();
                    for (index, &&(child, length)) in children.iter().enumerate().rev() {
                        tasks.push(Task::Length(length));
                        tasks.push(Task::Node(child, node));
                        if index > 0 {
                            tasks.push(Task::Text(","));
                        }
                    }
                }
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_quartets_and_branch_lengths() {
        let mut quartets = [0; 3];
        let mut length_sum = 0.0;
        for sample in 0..9000 {
            let tree = Tree::uniform(4, 0.1, sample_seed(1234, sample)).unwrap();
            assert_eq!(tree.edges.len(), 5);
            let mut neighbors = vec![Vec::new(); 6];
            for edge in &tree.edges {
                neighbors[edge.a].push(edge.b);
                neighbors[edge.b].push(edge.a);
                assert!(edge.length.is_finite() && edge.length > 0.0);
                length_sum += edge.length;
            }
            for &node in &neighbors[neighbors[3][0]] {
                if node < 3 {
                    quartets[node] += 1;
                }
            }
            for (node, adjacent) in neighbors.iter().enumerate() {
                assert_eq!(adjacent.len(), if node < 4 { 1 } else { 3 });
            }
        }
        assert!(quartets.iter().all(|&count| (2600..3400).contains(&count)));
        assert!((length_sum / 45000.0 - 0.1).abs() < 0.003);
    }

    #[test]
    fn sizes_and_validation() {
        for tips in [3, 20, 100, 10_000] {
            let tree = Tree::uniform(tips, 0.1, 42).unwrap();
            assert_eq!(tree.edges.len(), 2 * tips - 3);
            let text = tree.newick();
            assert_eq!(text.matches(':').count(), 2 * tips - 3);
            assert!(text.ends_with(";\n"));
        }
        assert!(Tree::uniform(2, 0.1, 1).is_err());
        assert!(Tree::uniform(4, f64::NAN, 1).is_err());
    }
}
