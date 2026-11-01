//! Application adapter around phylotree trees; sampling lives in that crate.
pub use crate::rng::sample_seed;
use phylotree::simulation::{self, BranchLengths, Topology};
use phylotree::tree::Tree as PhyloTree;
use rand::{SeedableRng, rngs::StdRng};

#[derive(Clone, Debug)]
pub struct Edge {
    pub a: usize,
    pub b: usize,
    pub length: f64,
}

#[derive(Debug)]
pub struct Tree {
    inner: PhyloTree,
    tips: usize,
    edges: Vec<Edge>,
}
impl Tree {
    pub fn tips(&self) -> usize {
        self.tips
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn phylogeny(&self) -> &PhyloTree {
        &self.inner
    }
    pub fn uniform(tips: usize, branch_mean: f64, seed: u64) -> Result<Self, String> {
        if !(3..=10_000).contains(&tips)
            || !branch_mean.is_finite()
            || !(1e-8..=100.0).contains(&branch_mean)
        {
            return Err("Need 3 <= tips <= 10000 and a finite branch mean in [1e-8, 100]".into());
        }
        let mut rng = StdRng::seed_from_u64(seed);
        Self::from_phylogeny(
            simulation::simulate(
                tips,
                Topology::UniformUnrooted,
                BranchLengths::Exponential { mean: branch_mean },
                &mut rng,
            )
            .map_err(|e| e.to_string())?,
        )
    }
    /// Map named tips and the display root to the plain edge ABI without Newick.
    pub fn from_phylogeny(inner: PhyloTree) -> Result<Self, String> {
        let root = inner.root().map_err(|e| e.to_string())?;
        let order = preorder(&inner)?;
        let tips = inner.n_leaves();
        if !(3..=10_000).contains(&tips) {
            return Err("AliSim adapter needs 3..10000 tips".into());
        }
        let mut map = std::collections::HashMap::new();
        let mut seen = vec![false; tips];
        for &id in &order {
            let node = inner.get(id).map_err(|e| e.to_string())?;
            if node.children.is_empty() {
                let label = node
                    .name
                    .as_deref()
                    .and_then(|n| n.strip_prefix('T'))
                    .and_then(|s| s.parse::<usize>().ok())
                    .filter(|&i| i > 0 && i <= tips)
                    .ok_or("Adapter expects unique labels T1..Tn")?
                    - 1;
                if seen[label] {
                    return Err("Duplicate tip label".into());
                }
                seen[label] = true;
                map.insert(id, label);
            }
        }
        map.insert(root, tips);
        let mut next = tips + 1;
        for &id in &order {
            if let std::collections::hash_map::Entry::Vacant(entry) = map.entry(id) {
                entry.insert(next);
                next += 1;
            }
        }
        let mut edges = Vec::new();
        for &id in &order {
            let node = inner.get(id).map_err(|e| e.to_string())?;
            if !(node.children.is_empty()
                || node.children.len() == 2
                || (id == root && node.children.len() == 3))
            {
                return Err("AliSim adapter expects a binary rooted or unrooted tree".into());
            }
            for &child in &node.children {
                let length = inner
                    .get(child)
                    .map_err(|e| e.to_string())?
                    .parent_edge
                    .ok_or("Missing branch length")?;
                if !length.is_finite() || length < 0.0 {
                    return Err("Invalid branch length".into());
                }
                edges.push(Edge {
                    a: map[&id],
                    b: map[&child],
                    length,
                });
            }
        }
        Ok(Self { inner, tips, edges })
    }
    pub fn newick(&self) -> String {
        format!(
            "{}\n",
            self.inner.to_newick().expect("validated generated tree")
        )
    }
}

/// Iterative preorder avoids recursive traversal on large/deep generated trees.
pub(crate) fn preorder(tree: &PhyloTree) -> Result<Vec<phylotree::tree::NodeId>, String> {
    let mut pending = vec![tree.root().map_err(|e| e.to_string())?];
    let mut order = Vec::new();
    while let Some(id) = pending.pop() {
        order.push(id);
        pending.extend(
            tree.get(id)
                .map_err(|e| e.to_string())?
                .children
                .iter()
                .rev()
                .copied(),
        );
    }
    Ok(order)
}
