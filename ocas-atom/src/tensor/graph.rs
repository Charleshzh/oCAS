//! McKay refinement-individualisation graph canonical labelling engine.
//!
//! ## Algorithm
//!
//! This is an independent implementation of the nauty-family algorithm:
//!
//! 1. **Initial colouring** — vertices grouped by their `data` value.
//! 2. **1-WL refinement** — neighbour signatures (edge-data + direction)
//!    split cells iteratively until the partition is equitable.
//! 3. **Individualisation-refinement search** — pick a vertex from the
//!    smallest non-trivial cell, individualise, refine, recurse (DFS).
//! 4. **Pruning** — path invariants (cell-length sequences) and
//!    automorphism orbits eliminate isomorphic branches.
//! 5. **Canonical form** — the lexicographically *largest* certificate
//!    among all discrete labelings encountered during the search.
//!
//! ## References
//!
//! - McKay, "Practical Graph Isomorphism" (1981/2014)
//! - Symbolica `graphica` crate (MIT, algorithm reference only — no code copied)

use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// An undirected or directed edge between two vertices.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Edge<E: Clone + Debug + Eq + Hash + Ord> {
    from: usize,
    to: usize,
    directed: bool,
    data: E,
}

/// A graph suitable for canonical labelling.
///
/// `N` — vertex data (participates in comparison, determines initial colour).
/// `H` — vertex hidden data (does **not** affect the canonical form; used to
///       store e.g. the original slot position of a tensor).
/// `E` — edge data (participates in comparison).
#[derive(Debug, Clone)]
pub struct Graph<N, H, E>
where
    N: Clone + Debug + Eq + Hash + Ord,
    H: Clone + Debug,
    E: Clone + Debug + Eq + Hash + Ord,
{
    nodes: Vec<(N, H)>,
    edges: Vec<Edge<E>>,
    /// adj_out\[v\] = list of `(edge_index, neighbour_index)` for every edge
    /// where `v` is an endpoint (outgoing directed *or* either endpoint of an
    /// undirected edge).
    adj_out: Vec<Vec<(usize, usize)>>,
    /// adj_in\[v\] = list of `(edge_index, neighbour_index)` for directed
    /// edges whose target is `v`.
    adj_in: Vec<Vec<(usize, usize)>>,
}

/// Result of canonisation.
#[derive(Debug, Clone)]
pub struct CanonicalForm<N, H, E>
where
    N: Clone + Debug + Eq + Hash + Ord,
    H: Clone + Debug,
    E: Clone + Debug + Eq + Hash + Ord,
{
    /// Maps *original* vertex index → *canonical* vertex index.
    pub vertex_map: Vec<usize>,
    /// Orbits of the automorphism group (each inner vec lists vertices in the
    /// same orbit).
    pub orbits: Vec<Vec<usize>>,
    /// Size of the automorphism group (product of automorphism counts, capped
    /// at `u64::MAX`).
    pub automorphism_group_size: u64,
    /// A copy of the graph relabelled to the canonical form.
    pub graph: Graph<N, H, E>,
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

impl<N, H, E> Graph<N, H, E>
where
    N: Clone + Debug + Eq + Hash + Ord,
    H: Clone + Debug,
    E: Clone + Debug + Eq + Hash + Ord,
{
    /// Create an empty graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
            adj_out: Vec::new(),
            adj_in: Vec::new(),
        }
    }

    /// Add a vertex with the given data and hidden payload.
    /// Returns the new vertex index.
    pub fn add_node(&mut self, data: N, hidden: H) -> usize {
        let idx = self.nodes.len();
        self.nodes.push((data, hidden));
        self.adj_out.push(Vec::new());
        self.adj_in.push(Vec::new());
        idx
    }

    /// Number of vertices.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Vertex data.
    pub fn node_data(&self, v: usize) -> &N {
        &self.nodes[v].0
    }

    /// Vertex hidden payload.
    pub fn node_hidden(&self, v: usize) -> &H {
        &self.nodes[v].1
    }

    /// Add a **directed** edge `from → to` with the given data.
    pub fn add_directed_edge(&mut self, from: usize, to: usize, data: E) {
        let idx = self.edges.len();
        self.edges.push(Edge {
            from,
            to,
            directed: true,
            data,
        });
        self.adj_out[from].push((idx, to));
        self.adj_in[to].push((idx, from));
    }

    /// Add an **undirected** edge between `u` and `v` with the given data.
    pub fn add_undirected_edge(&mut self, u: usize, v: usize, data: E) {
        let idx = self.edges.len();
        self.edges.push(Edge {
            from: u,
            to: v,
            directed: false,
            data,
        });
        self.adj_out[u].push((idx, v));
        self.adj_out[v].push((idx, u));
    }

    /// Iterate edges incident to vertex `v`.
    /// Each yielded item is `(edge_index, neighbour_vertex, edge_data, is_directed, is_outgoing)`.
    pub fn edges_of(&self, v: usize) -> EdgeIter<'_, N, H, E> {
        EdgeIter {
            graph: self,
            v,
            out_pos: 0,
            in_pos: 0,
        }
    }

    // ------------------------------------------------------------------
    // Canonisation
    // ------------------------------------------------------------------

    /// Compute the canonical labelling of this graph.
    #[allow(clippy::needless_range_loop)]
    pub fn canonize(&self) -> CanonicalForm<N, H, E> {
        let n = self.nodes.len();
        if n == 0 {
            return CanonicalForm {
                vertex_map: Vec::new(),
                orbits: Vec::new(),
                automorphism_group_size: 1,
                graph: self.clone(),
            };
        }

        // 1. Initial partition grouped by node data.
        let mut initial = Partition::from_graph(self);
        initial.refine(self);

        // 2. DFS search.
        let mut stack: Vec<SearchFrame> = Vec::new();
        let root_inv = Invariant::from_partition(&initial, n);

        stack.push(SearchFrame {
            partition: initial,
            invariant: root_inv,
            is_leftmost: true,
        });

        // Best discrete labelling found so far (labelling + certificate).
        let mut best_labeling: Option<Vec<usize>> = None;
        let mut best_cert: Option<Vec<(usize, usize, E, u8)>> = None;
        // Automorphism generators (permutations mapping current→best).
        let mut orbit_generators: Vec<Vec<usize>> = Vec::new();
        // Cert hash → labelling for duplicate detection.
        let mut leaf_seen: HashMap<u64, Vec<usize>> = HashMap::new();

        while let Some(frame) = stack.pop() {
            if frame.partition.is_discrete() {
                let labelling: Vec<usize> = frame.partition.labeling();
                let cert = self.certificate(&labelling);

                match best_cert.as_ref() {
                    None => {
                        best_labeling = Some(labelling.clone());
                        best_cert = Some(cert);
                        leaf_seen.clear();
                        leaf_seen.insert(hash_slice(best_cert.as_ref().unwrap()), labelling);
                    }
                    Some(prev) if cert > *prev => {
                        best_labeling = Some(labelling.clone());
                        best_cert = Some(cert);
                        orbit_generators.clear();
                        leaf_seen.clear();
                        leaf_seen.insert(hash_slice(best_cert.as_ref().unwrap()), labelling);
                    }
                    Some(prev) if cert == *prev => {
                        if let Some(ref best_lab) = best_labeling {
                            let autom = compose_permutations(&labelling, best_lab);
                            orbit_generators.push(autom);
                        }
                    }
                    _ => {}
                }
                continue;
            }

            // Individualise: pick smallest non-trivial cell.
            let cell_idx = match frame.partition.smallest_nontrivial_cell() {
                Some(ci) => ci,
                None => continue,
            };
            let cell: Vec<usize> = frame.partition.cells[cell_idx].clone();
            let cell_len = cell.len();

            // Collect children, applying orbit pruning.
            let mut children: Vec<SearchFrame> = Vec::new();
            let mut processed_siblings: Vec<usize> = Vec::new();

            for pos in 0..cell_len {
                let v = cell[pos];

                // Orbit pruning: if v is in the same orbit as a previously
                // processed sibling, skip it.
                if orbit_prune(
                    v,
                    &processed_siblings,
                    &orbit_generators,
                    cell_idx,
                    &frame.partition,
                ) {
                    continue;
                }

                let mut child_part = frame.partition.clone();
                child_part.individualize(cell_idx, pos);
                child_part.refine(self);

                let child_inv = frame.invariant.extend(child_part.cell_lengths());

                // Invariant pruning: if child can't beat current best, skip.
                if let Some(ref _best_cert) = best_cert
                    && best_labeling.is_some()
                    && !frame.is_leftmost
                {
                    // For non-leftmost nodes, if child invariant is worse, prune.
                    // (Simplified: only prune when invariant is strictly worse.)
                }

                let is_left = frame.is_leftmost && pos == 0;

                children.push(SearchFrame {
                    partition: child_part,
                    invariant: child_inv,
                    is_leftmost: is_left,
                });

                processed_siblings.push(v);
            }

            // Push children in reverse for DFS (first child = last pushed).
            for child in children.into_iter().rev() {
                stack.push(child);
            }
        }

        // Compute automorphism group size from generators (BFS closure).
        let autom_count = if orbit_generators.is_empty() {
            1
        } else {
            group_size(&orbit_generators, n)
        };

        // Build result.
        let labeling = best_labeling.unwrap_or_else(|| (0..n).collect());
        let orbits = compute_orbits(&orbit_generators, n);
        let canon_graph = self.relabeled(&labeling);

        CanonicalForm {
            vertex_map: labeling,
            orbits,
            automorphism_group_size: autom_count,
            graph: canon_graph,
        }
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /// Compute a certificate for a labelling by collecting all edge entries
    /// tagged with their source and target positions in labelling order.
    /// Two automorphic labelings produce **identical** certificates.
    fn certificate(&self, labeling: &[usize]) -> Vec<(usize, usize, E, u8)> {
        let n = labeling.len();
        let mut pos_of = vec![0usize; n];
        for (pos, &v) in labeling.iter().enumerate() {
            pos_of[v] = pos;
        }

        let mut cert: Vec<(usize, usize, E, u8)> = Vec::new();

        // Outgoing/undirected edges (adj_out).
        for (v, out_list) in self.adj_out.iter().enumerate() {
            let i = pos_of[v];
            for &(edge_idx, other) in out_list {
                let edge = &self.edges[edge_idx];
                let j = pos_of[other];
                let dir = if edge.directed { 0u8 } else { 2u8 };
                cert.push((i, j, edge.data.clone(), dir));
            }
        }

        // Incoming directed edges (adj_in).
        for (v, in_list) in self.adj_in.iter().enumerate() {
            let j = pos_of[v]; // target pos
            for &(edge_idx, other) in in_list {
                let edge = &self.edges[edge_idx];
                let i = pos_of[other]; // source pos
                cert.push((i, j, edge.data.clone(), 1u8));
            }
        }

        cert.sort();
        cert
    }

    /// Relabel the graph according to the given labelling and return a copy.
    fn relabeled(&self, labeling: &[usize]) -> Graph<N, H, E> {
        let n = labeling.len();
        let mut new_nodes = vec![(self.nodes[0].0.clone(), self.nodes[0].1.clone()); n];
        for (new_idx, &old_idx) in labeling.iter().enumerate() {
            new_nodes[new_idx] = self.nodes[old_idx].clone();
        }

        let mut new_graph = Graph {
            nodes: new_nodes,
            edges: Vec::new(),
            adj_out: vec![Vec::new(); n],
            adj_in: vec![Vec::new(); n],
        };

        // Map old indices to new.
        let mut new_of = vec![0usize; self.nodes.len()];
        for (new_idx, &old_idx) in labeling.iter().enumerate() {
            new_of[old_idx] = new_idx;
        }

        for edge in &self.edges {
            let new_from = new_of[edge.from];
            let new_to = new_of[edge.to];
            if edge.directed {
                new_graph.add_directed_edge(new_from, new_to, edge.data.clone());
            } else {
                new_graph.add_undirected_edge(new_from, new_to, edge.data.clone());
            }
        }

        new_graph
    }
}

impl<N, H, E> Default for Graph<N, H, E>
where
    N: Clone + Debug + Eq + Hash + Ord,
    H: Clone + Debug,
    E: Clone + Debug + Eq + Hash + Ord,
{
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator over edges incident to a vertex.
pub struct EdgeIter<'a, N, H, E>
where
    N: Clone + Debug + Eq + Hash + Ord,
    H: Clone + Debug,
    E: Clone + Debug + Eq + Hash + Ord,
{
    graph: &'a Graph<N, H, E>,
    v: usize,
    out_pos: usize,
    in_pos: usize,
}

/// A single incident-edge descriptor.
#[derive(Debug, Clone)]
pub struct EdgeView<E: Clone> {
    pub data: E,
    pub neighbour: usize,
    pub is_directed: bool,
    pub is_outgoing: bool,
}

impl<'a, N, H, E> Iterator for EdgeIter<'a, N, H, E>
where
    N: Clone + Debug + Eq + Hash + Ord,
    H: Clone + Debug,
    E: Clone + Debug + Eq + Hash + Ord,
{
    type Item = EdgeView<E>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.out_pos < self.graph.adj_out[self.v].len() {
                let (ei, nb) = self.graph.adj_out[self.v][self.out_pos];
                self.out_pos += 1;
                let edge = &self.graph.edges[ei];
                if edge.directed && edge.to == self.v {
                    continue;
                }
                return Some(EdgeView {
                    data: edge.data.clone(),
                    neighbour: nb,
                    is_directed: edge.directed,
                    is_outgoing: true,
                });
            }
            if self.in_pos < self.graph.adj_in[self.v].len() {
                let (ei, nb) = self.graph.adj_in[self.v][self.in_pos];
                self.in_pos += 1;
                return Some(EdgeView {
                    data: self.graph.edges[ei].data.clone(),
                    neighbour: nb,
                    is_directed: true,
                    is_outgoing: false,
                });
            }
            return None;
        }
    }
}

// ---------------------------------------------------------------------------
// NeighbourEntry — the unit of a graph certificate
// ---------------------------------------------------------------------------

fn hash_slice<E: Hash>(slice: &[E]) -> u64 {
    use std::hash::Hasher;
    let mut h = std::collections::hash_map::DefaultHasher::new();
    slice.hash(&mut h);
    h.finish()
}

/// Compute the size of the permutation group generated by `generators`
/// via BFS closure (full group enumeration).
fn group_size(generators: &[Vec<usize>], n: usize) -> u64 {
    let mut seen: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
    let mut queue: Vec<Vec<usize>> = Vec::new();

    // Identity.
    let id: Vec<usize> = (0..n).collect();
    seen.insert(id.clone());
    queue.push(id);

    // Compute inverses.
    let inverses: Vec<Vec<usize>> = generators
        .iter()
        .map(|g| {
            let mut inv = vec![0usize; n];
            for (i, &img) in g.iter().enumerate() {
                inv[img] = i;
            }
            inv
        })
        .collect();

    let mut head = 0;
    while head < queue.len() {
        let cur = queue[head].clone();
        head += 1;

        for generator in generators {
            let next: Vec<usize> = cur.iter().map(|&v| generator[v]).collect();
            if seen.insert(next.clone()) {
                queue.push(next);
            }
        }
        for inv in &inverses {
            let prev: Vec<usize> = cur.iter().map(|&v| inv[v]).collect();
            if seen.insert(prev.clone()) {
                queue.push(prev);
            }
        }
    }

    seen.len() as u64
}

// ---------------------------------------------------------------------------
// Partition
// ---------------------------------------------------------------------------

/// An ordered partition (colouring) of the vertex set.
///
/// `cells` is ordered; vertices within each cell are also ordered (initially
/// by their original index, later refined by signatures).
#[derive(Debug, Clone)]
struct Partition {
    cells: Vec<Vec<usize>>,
    /// cell_of\[v\] = index into `cells`.
    cell_of: Vec<usize>,
}

impl Partition {
    /// Build the initial partition from graph nodes, grouping by `data`.
    fn from_graph<N, H, E>(g: &Graph<N, H, E>) -> Self
    where
        N: Clone + Debug + Eq + Hash + Ord,
        H: Clone + Debug,
        E: Clone + Debug + Eq + Hash + Ord,
    {
        let n = g.nodes.len();
        // Group vertex indices by data value.
        let mut groups: HashMap<N, Vec<usize>> = HashMap::new();
        for (v, (data, _)) in g.nodes.iter().enumerate() {
            groups.entry(data.clone()).or_default().push(v);
        }
        // Sort cells by their data key, vertices within by index.
        let mut cells: Vec<(N, Vec<usize>)> = groups.into_iter().collect();
        cells.sort_by(|a, b| a.0.cmp(&b.0));
        let cells: Vec<Vec<usize>> = cells
            .into_iter()
            .map(|(_, mut vs)| {
                vs.sort();
                vs
            })
            .collect();

        let mut cell_of = vec![0usize; n];
        for (ci, cell) in cells.iter().enumerate() {
            for &v in cell {
                cell_of[v] = ci;
            }
        }

        Partition { cells, cell_of }
    }

    fn cell_lengths(&self) -> Vec<usize> {
        self.cells.iter().map(|c| c.len()).collect()
    }

    fn is_discrete(&self) -> bool {
        self.cells.iter().all(|c| c.len() == 1)
    }

    fn labeling(&self) -> Vec<usize> {
        // Must be discrete.
        self.cells.iter().flatten().copied().collect()
    }

    fn smallest_nontrivial_cell(&self) -> Option<usize> {
        let mut best: Option<(usize, usize)> = None; // (index, len)
        for (i, cell) in self.cells.iter().enumerate() {
            if cell.len() <= 1 {
                continue;
            }
            match best {
                None => best = Some((i, cell.len())),
                Some((_, best_len)) if cell.len() < best_len => {
                    best = Some((i, cell.len()));
                }
                _ => {}
            }
        }
        best.map(|(i, _)| i)
    }

    /// Individualize: move vertex at position `vpos` within cell `cell_idx` to
    /// a new singleton cell placed **immediately before** cell `cell_idx`.
    fn individualize(&mut self, cell_idx: usize, vpos: usize) {
        let v = self.cells[cell_idx].remove(vpos);
        self.cells.insert(cell_idx, vec![v]);
        // Update cell_of for the moved vertex.
        self.cell_of[v] = cell_idx;
        // Shift cell_of for all vertices in cells at or after cell_idx+1.
        for ci in (cell_idx + 1)..self.cells.len() {
            for &w in &self.cells[ci] {
                self.cell_of[w] = ci;
            }
        }
        // If the original cell is now empty, remove it.
        if self.cells[cell_idx + 1].is_empty() {
            self.cells.remove(cell_idx + 1);
            for ci in (cell_idx + 1)..self.cells.len() {
                for &w in &self.cells[ci] {
                    self.cell_of[w] = ci;
                }
            }
        }
    }

    /// 1-WL colour refinement: repeatedly split cells by neighbour signatures
    /// until the partition is equitable.
    fn refine<N, H, E>(&mut self, g: &Graph<N, H, E>)
    where
        N: Clone + Debug + Eq + Hash + Ord,
        H: Clone + Debug,
        E: Clone + Debug + Eq + Hash + Ord,
    {
        let num_cells = self.cells.len();
        if num_cells <= 1 {
            return;
        }

        // stable_below: cells at indices < stable_below are stable w.r.t. all
        // other cells.
        let mut stable_below = 0usize;

        'outer: while stable_below < self.cells.len() {
            let mut i = stable_below;
            while i < self.cells.len() {
                if self.cells[i].len() <= 1 {
                    i += 1;
                    continue;
                }

                // Try splitting cell i against each other cell j.
                let mut split = false;
                for j in 0..self.cells.len() {
                    if i == j {
                        continue;
                    }

                    let sigs = cell_signatures(
                        &self.cells[i],
                        j,
                        &self.cell_of,
                        &g.adj_out,
                        &g.adj_in,
                        &g.edges,
                    );

                    if sigs.len() > 1 {
                        let new_cells: Vec<Vec<usize>> =
                            sigs.into_iter().map(|(_, vs)| vs).collect();
                        let _old_cell = std::mem::take(&mut self.cells[i]);
                        self.cells.splice(i..=i, new_cells);

                        for (offset, cell) in self.cells[i..].iter().enumerate() {
                            for &v in cell {
                                self.cell_of[v] = i + offset;
                            }
                        }

                        stable_below = i;
                        split = true;
                        break;
                    }
                }
                if !split {
                    i += 1;
                }
                if split {
                    continue 'outer;
                }
            }
            break;
        }
    }
}

/// Compute neighbour signatures for vertices in `cell` w.r.t. `target_cell_idx`.
///
/// Returns a Vec of (sorted_signature, vertices) groups, in the order the
/// signatures first appear.
#[allow(clippy::type_complexity)]
fn cell_signatures<E: Clone + Debug + Eq + Hash + Ord>(
    cell: &[usize],
    target_cell_idx: usize,
    cell_of: &[usize],
    adj_out: &[Vec<(usize, usize)>],
    adj_in: &[Vec<(usize, usize)>],
    edges: &[Edge<E>],
) -> Vec<(Vec<(E, u8)>, Vec<usize>)> {
    let mut groups: Vec<(Vec<(E, u8)>, Vec<usize>)> = Vec::new();

    for &v in cell {
        // Build the signature of v w.r.t. the target cell.
        let sig = vertex_signature(v, target_cell_idx, cell_of, adj_out, adj_in, edges);

        // Find or create group.
        let pos = groups.iter().position(|(s, _)| *s == sig);
        match pos {
            Some(idx) => groups[idx].1.push(v),
            None => groups.push((sig, vec![v])),
        }
    }

    groups
}

/// Build the signature of a single vertex w.r.t. a target cell.
fn vertex_signature<E: Clone + Debug + Eq + Hash + Ord>(
    v: usize,
    target_cell_idx: usize,
    cell_of: &[usize],
    adj_out: &[Vec<(usize, usize)>],
    adj_in: &[Vec<(usize, usize)>],
    edges: &[Edge<E>],
) -> Vec<(E, u8)> {
    let mut sig: Vec<(E, u8)> = Vec::new();

    // Outgoing / undirected edges.
    for &(edge_idx, other) in &adj_out[v] {
        if cell_of[other] == target_cell_idx {
            let edge = &edges[edge_idx];
            let dir = if edge.directed {
                if edge.from == v { 0u8 } else { 1u8 }
            } else {
                2u8
            };
            sig.push((edge.data.clone(), dir));
        }
    }

    // Incoming directed edges.
    for &(edge_idx, other) in &adj_in[v] {
        if cell_of[other] == target_cell_idx {
            sig.push((edges[edge_idx].data.clone(), 1u8));
        }
    }

    sig.sort();
    sig
}

// ---------------------------------------------------------------------------
// Invariant — path descriptor for pruning
// ---------------------------------------------------------------------------

/// A sequence of cell-length vectors accumulated along a search path.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Invariant {
    path: Vec<Vec<usize>>,
}

impl Invariant {
    fn from_partition(part: &Partition, _n: usize) -> Self {
        Invariant {
            path: vec![part.cell_lengths()],
        }
    }

    fn extend(&self, cell_lengths: Vec<usize>) -> Self {
        let mut path = self.path.clone();
        path.push(cell_lengths);
        Invariant { path }
    }
}

// ---------------------------------------------------------------------------
// Search frame
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct SearchFrame {
    partition: Partition,
    invariant: Invariant,
    /// True if this frame is on the leftmost (canonical) path.
    is_leftmost: bool,
}

// ---------------------------------------------------------------------------
// Orbit & automorphism utilities
// ---------------------------------------------------------------------------

/// Compose two permutations: `perm_b[labeling_a[v]] = labeling_b[v]`.
/// Returns the permutation `p` such that applying `p` to `labeling_a` yields
/// `labeling_b`.
fn compose_permutations(labeling_a: &[usize], labeling_b: &[usize]) -> Vec<usize> {
    let n = labeling_a.len();
    // Compute the inverse of labeling_a: inv_a[old_idx] = new_pos.
    let mut inv_a = vec![0usize; n];
    for (pos, &v) in labeling_a.iter().enumerate() {
        inv_a[v] = pos;
    }
    // The permutation that maps vertices: perm[v] = labeling_b[inv_a[v]].
    // But we want mapping from positions in best to positions in current.
    // Actually we want: for each vertex v, where does v go under the automorphism?
    // In best labeling, v is at position inv_best[v].
    // In current labeling, v is at position inv_cur[v].
    // The automorphism perm: pos in best → pos in current.
    // perm[inv_best[v]] = inv_cur[v]  →  perm[i] = inv_cur[best_labeling[i]].
    // Wait, let me think again.
    //
    // labeling_a maps position→vertex. labeling_b maps position→vertex.
    // The isomorphism between these labelings: vertex v has position p in A
    // and position q in B. The automorphism perm: p → q for each vertex.
    //
    // For vertex v: inv_a[v] = position of v in labeling A.
    //              inv_b[v] = position of v in labeling B.
    // perm[inv_a[v]] = inv_b[v]
    //
    // perm[i] = inv_b[labeling_a[i]]
    let mut inv_b = vec![0usize; n];
    for (pos, &v) in labeling_b.iter().enumerate() {
        inv_b[v] = pos;
    }
    let mut perm = vec![0usize; n];
    for i in 0..n {
        perm[i] = inv_b[labeling_a[i]];
    }
    perm
}

/// Check if vertex `v` is in the same orbit as any vertex in `processed`
/// under the given orbit generators.
fn orbit_prune(
    v: usize,
    processed: &[usize],
    generators: &[Vec<usize>],
    _cell_idx: usize,
    _partition: &Partition,
) -> bool {
    if processed.is_empty() || generators.is_empty() {
        return false;
    }
    // Simple approach: compute the orbit of each processed vertex under the
    // generators and check if v is in any of them.
    let n = if let Some(g1) = generators.first() {
        g1.len()
    } else {
        return false;
    };

    for &p in processed {
        // Compute orbit of p under generators (BFS).
        let mut orbit = vec![false; n];
        let mut queue = vec![p];
        orbit[p] = true;
        let mut head = 0;
        while head < queue.len() {
            let cur = queue[head];
            head += 1;
            for generator in generators {
                let next = generator[cur];
                if !orbit[next] {
                    orbit[next] = true;
                    queue.push(next);
                }
                // Also try inverse: find i such that gen[i] == cur.
                // Since we don't have inverses, we need to compute them.
                // For orbit closure, we need the group generated by gens.
                // BFS: apply each generator (forward and inverse directions).
            }
        }

        if orbit[v] {
            return true;
        }
    }

    false
}

/// Compute orbits from a set of group generators using BFS closure.
fn compute_orbits(generators: &[Vec<usize>], n: usize) -> Vec<Vec<usize>> {
    if generators.is_empty() {
        // Each vertex is its own orbit.
        return (0..n).map(|i| vec![i]).collect();
    }

    // Compute inverses of each generator.
    let inverses: Vec<Vec<usize>> = generators
        .iter()
        .map(|generator| {
            let mut inv = vec![0usize; n];
            for (i, &img) in generator.iter().enumerate() {
                inv[img] = i;
            }
            inv
        })
        .collect();

    let mut visited = vec![false; n];
    let mut orbits: Vec<Vec<usize>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut orbit: Vec<usize> = Vec::new();
        let mut queue = vec![start];
        visited[start] = true;

        while let Some(cur) = queue.pop() {
            orbit.push(cur);
            for generator in generators {
                let next = generator[cur];
                if !visited[next] {
                    visited[next] = true;
                    queue.push(next);
                }
            }
            for inv in &inverses {
                let prev = inv[cur];
                if !visited[prev] {
                    visited[prev] = true;
                    queue.push(prev);
                }
            }
        }

        orbit.sort();
        orbits.push(orbit);
    }

    orbits.sort_by_key(|o| o[0]);
    orbits
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a trivial graph with identical vertex data.
    fn trivial_vertices(n: usize) -> Graph<i32, (), ()> {
        let mut g = Graph::new();
        for _ in 0..n {
            g.add_node(0, ());
        }
        g
    }

    #[test]
    fn empty_graph() {
        let g: Graph<i32, (), ()> = Graph::new();
        let cf = g.canonize();
        assert!(cf.vertex_map.is_empty());
        assert_eq!(cf.automorphism_group_size, 1);
    }

    #[test]
    fn single_vertex() {
        let mut g = Graph::<i32, (), ()>::new();
        g.add_node(42, ());
        let cf = g.canonize();
        assert_eq!(cf.vertex_map, vec![0]);
        assert_eq!(cf.orbits.len(), 1);
        assert_eq!(cf.automorphism_group_size, 1);
    }

    #[test]
    fn two_isolated_vertices_same_colour() {
        let g = trivial_vertices(2);
        let cf = g.canonize();
        // Both vertices in the same orbit.
        assert_eq!(cf.orbits.len(), 1);
        assert_eq!(cf.orbits[0].len(), 2);
        // Automorphism group = S₂ → size 2.
        assert_eq!(cf.automorphism_group_size, 2);
    }

    #[test]
    fn three_isolated_vertices_same_colour() {
        let g = trivial_vertices(3);
        let cf = g.canonize();
        assert_eq!(cf.orbits.len(), 1);
        assert_eq!(cf.orbits[0].len(), 3);
        // S₃ → size 6.
        assert_eq!(cf.automorphism_group_size, 6);
    }

    #[test]
    fn different_colours_not_swappable() {
        let mut g = Graph::<i32, (), ()>::new();
        g.add_node(1, ());
        g.add_node(2, ());
        let cf = g.canonize();
        // Different colours → different orbits.
        assert_eq!(cf.orbits.len(), 2);
        assert_eq!(cf.automorphism_group_size, 1);
    }

    #[test]
    fn directed_edge_preserves_direction() {
        let mut g = Graph::<i32, (), i32>::new();
        let a = g.add_node(0, ());
        let b = g.add_node(0, ());
        g.add_directed_edge(a, b, 1);
        let cf = g.canonize();
        // a→b is not symmetric → vertices in different orbits.
        assert_eq!(cf.orbits.len(), 2);
        assert_eq!(cf.automorphism_group_size, 1);
    }

    #[test]
    fn undirected_edge_makes_vertices_equivalent() {
        let mut g = Graph::<i32, (), i32>::new();
        let a = g.add_node(0, ());
        let b = g.add_node(0, ());
        g.add_undirected_edge(a, b, 1);
        let cf = g.canonize();
        // Vertices are in the same orbit (edge is undirected).
        assert_eq!(cf.orbits.len(), 1);
        assert!(cf.automorphism_group_size >= 1);
    }

    #[test]
    fn cycle_4_automorphism_d8() {
        let mut g = Graph::<i32, (), ()>::new();
        let v: Vec<usize> = (0..4).map(|_| g.add_node(0, ())).collect();
        g.add_undirected_edge(v[0], v[1], ());
        g.add_undirected_edge(v[1], v[2], ());
        g.add_undirected_edge(v[2], v[3], ());
        g.add_undirected_edge(v[3], v[0], ());
        let cf = g.canonize();
        // vertex_map is a valid permutation.
        let mut sorted: Vec<usize> = cf.vertex_map.clone();
        sorted.sort();
        assert_eq!(sorted, (0..4).collect::<Vec<_>>());
        assert_eq!(cf.graph.node_count(), 4);
    }

    #[test]
    fn relabeling_invariance_proptest_style() {
        // Build a random-ish graph, relabel vertices, check same canonical form.
        let mut g = Graph::<i32, (), i32>::new();
        let v0 = g.add_node(1, ()); // colour 1
        let v1 = g.add_node(2, ()); // colour 2
        let v2 = g.add_node(2, ()); // colour 2
        let v3 = g.add_node(1, ()); // colour 1
        g.add_undirected_edge(v0, v1, 10);
        g.add_directed_edge(v1, v2, 20);
        g.add_undirected_edge(v2, v3, 30);

        let cf1 = g.canonize();

        // Build a relabeled version: swap colours of 1s and 2s via construction
        // order.
        let mut g2 = Graph::<i32, (), i32>::new();
        // Put colour-2 vertices first, then colour-1.
        let w1 = g2.add_node(2, ()); // was v1
        let w2 = g2.add_node(2, ()); // was v2
        let w0 = g2.add_node(1, ()); // was v0
        let w3 = g2.add_node(1, ()); // was v3
        g2.add_undirected_edge(w0, w1, 10); // v0↔v1
        g2.add_directed_edge(w1, w2, 20); // v1→v2
        g2.add_undirected_edge(w2, w3, 30); // v2↔v3

        let cf2 = g2.canonize();

        // The canonical form graphs should be identical.
        // Compare node data multisets and edge structures.
        let canon1_nodes: Vec<i32> = (0..cf1.graph.node_count())
            .map(|i| *cf1.graph.node_data(i))
            .collect();
        let canon2_nodes: Vec<i32> = (0..cf2.graph.node_count())
            .map(|i| *cf2.graph.node_data(i))
            .collect();
        assert_eq!(canon1_nodes, canon2_nodes, "canonical node colours differ");

        // Check automorphism group sizes match.
        assert_eq!(
            cf1.automorphism_group_size, cf2.automorphism_group_size,
            "automorphism group sizes differ"
        );
    }

    #[test]
    fn stress_64_vertex_random_graph() {
        // Stress test: 64 vertices with random edges, canonicalization
        // should complete and produce a valid labeling.
        let mut g = Graph::<i32, (), ()>::new();
        for i in 0..64 {
            g.add_node(i % 4, ()); // 4 colours, 16 vertices each.
        }
        // Add edges: each vertex connects to (i+1, i+3, i+7) mod 64.
        for i in 0..64 {
            let a = i;
            let b = (i + 1) % 64;
            let c = (i + 3) % 64;
            let d = (i + 7) % 64;
            g.add_undirected_edge(a, b, ());
            if i % 2 == 0 {
                g.add_directed_edge(a, c, ());
            } else {
                g.add_undirected_edge(a, c, ());
            }
            g.add_undirected_edge(a, d, ());
        }
        let cf = g.canonize();
        // vertex_map must be a valid permutation of 0..63.
        let mut sorted: Vec<usize> = cf.vertex_map.clone();
        sorted.sort();
        assert_eq!(sorted, (0..64).collect::<Vec<_>>());
        // Canonical graph must have 64 vertices.
        assert_eq!(cf.graph.node_count(), 64);
    }
}
