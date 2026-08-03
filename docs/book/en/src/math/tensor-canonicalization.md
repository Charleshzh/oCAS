# Tensor Algebra & Canonicalization

This chapter discusses the algebraic structure of tensor expressions in oCAS and the mathematical foundations of canonicalization. The central question: given a tensor-product expression with contractions and symmetries, how can one decide whether two seemingly different writings (for example $T_{\mu\nu} U^{\nu\rho}$ and $U^{\rho\nu} T_{\nu\mu}$) represent the same mathematical object? oCAS's approach is to encode the tensor expression as a colored directed graph and then compute its **canonical form** via a graph-isomorphism algorithm.

---

## Prerequisites

### Tensor Products and Contractions

A tensor $T$ of type $(r, s)$ has $r$ **contravariant indices** (superscripts) and $s$ **covariant indices** (subscripts). The **tensor product** $T \otimes U$ of two tensors simply concatenates their indices.

**Contraction** pairs one contravariant index with one covariant index and sums over it. Under the Einstein summation convention:

$$
T_{\mu\nu} U^{\nu\rho} = \sum_{\nu} T_{\mu\nu} U^{\nu\rho}
$$

Here $\nu$ occurs twice — once as a subscript and once as a superscript — and is therefore summed away; the result is a tensor of type $(1,1)$. The summed-away indices are called **dummy indices** and the remaining ones **free indices**.

> **Einstein summation vs. explicit matching in oCAS**: the Einstein convention sums implicitly over repeated indices and relies on positions (upper/lower) to determine pairings. oCAS uses explicit matching — every index slot is an argument carrying a position marker, and a contraction is expressed by two slots sharing the same label. This removes the implicit dependence on "upper/lower" positions and makes the system more explicit.

### Raising and Lowering Indices

In a space with a metric $g_{\mu\nu}$, indices can be **raised** and **lowered** by contracting with the metric:

$$
T^{\mu} = g^{\mu\nu} T_{\nu}, \qquad T_{\mu} = g_{\mu\nu} T^{\nu}
$$

In the oCAS tensor system, index raising/lowering is marked on each index slot by the `IndexPosition` enum (`Upper` / `Lower`); no metric tensor is introduced automatically.

---

## Basic Concepts

### Index Matching

Each argument of a tensor function $T(i, j, k)$ is called an **index slot**. Every slot carries:

- a **label**: an `Atom`, such as `i`, `mu`, `d0`
- a **variance**: `Upper` (contravariant) or `Lower` (covariant)

A **contraction** between two slots is expressed by label matching: when two slots share the same label and one is `Upper` while the other is `Lower`, they form a contraction pair.

### Dummy-Index Management

In a tensor product $T(i, j) \cdot U(j, k)$, the label $j$ occurs exactly twice and is therefore a dummy index. oCAS's dummy-index management follows these rules:

1. **Identification**: walk the arguments of all `Fun` nodes in the expression and count the occurrences of every label. Labels occurring exactly **2** times are dummy indices.
2. **Validation**: in canonicalization (`canon.rs`), a label occurring more than twice raises `TensorCanonError::ContractedMoreThanOnce`; `dummy.rs` declares the `OverContracted` and `BadContraction` (same-variance pair) error variants, but `refresh_dummies` itself does not perform a variance check.
3. **Renaming**: after canonicalization, dummy indices are replaced by new names `d0`, `d1`, `d2`, … numbered per group. Different index groups (e.g. spacetime indices and internal indices) use different namespaces `d{group}_{n}`.

```rust
// Example of dummy-index renaming
// Input: T(i, j) * U(j, k) * V(i, l)
// Output: T(d0, d1) * U(d1, k) * V(d0, l)
// i → d0 (group 0, 0th), j → d1 (group 0, 1st)
```

This renaming is implemented by the `refresh_dummies` function in `dummy.rs`. It queries the index group to which a label belongs via `TensorRegistry::index_group`, ensuring that dummy indices of different dimensions use different namespaces.

---

## Core Theory

### Graph Encoding

The first step of tensor canonicalization is to encode the tensor expression as a colored directed graph $G = (V, E)$. The encoding rules are as follows (implemented by `tensor_to_graph` in `canon.rs`):

#### Vertex Types

| Type | Vertex data `TgNode` | Meaning |
|---|---|---|
| Head vertex | `Head(hash)` | a tensor function head $T$; `hash` is the hash of the function name |
| Slot vertex | `Slot(hash)` | an index slot; `hash` is the hash of the index label (symmetric slots uniformly use `Slot(0)`) |
| Scalar vertex | `Scalar(hash)` | a scalar factor in the product |

#### Edge Types

| Type | Edge data `TgEdge` | Direction | Meaning |
|---|---|---|---|
| Head→slot | `HeadToSlot(pos, flag)` | directed | from a tensor head to its `pos`-th slot; `flag=1` visible, `flag=0` hidden (symmetric) |
| Contraction | `Contraction(group)` | undirected | an edge between two contracted slot vertices; `group` identifies the index group |

#### Handling Symmetries

When some slots of a tensor are symmetric (`SymmetrySpec::symmetric_subsets`), these slots are marked **hidden** in the graph:

- The vertex data of symmetric slots is uniformly `Slot(0)`, so the graph-isomorphism engine may freely permute them.
- The head→slot edges of symmetric slots are marked `HeadToSlot(pos, 0)` (flag=0) and use the **sorted position** instead of the original one, so that $T(a, b)$ and $T(b, a)$ produce identical encodings (the flag itself still participates in the certificate comparison).
- At input time, symmetric slots are pre-sorted by the alphabetical order of their labels, ensuring that $T(a, b)$ and $T(b, a)$ generate the same graph.

**Key invariant**: two tensor expressions are mathematically equivalent $\Leftrightarrow$ their graphs are isomorphic.

### The McKay Refinement–Individualization Algorithm

Graph-isomorphism testing and canonical-labeling computation use the McKay refinement–individualization algorithm (of the nauty family), implemented in `graph.rs` (a standalone implementation of about 1100 lines).

#### Algorithm Overview

Given a graph $G$, the algorithm outputs:

1. **Canonical labeling** $\lambda: V \to \{0, 1, \dots, n-1\}$ — a bijection such that for any isomorphic graph $G'$, the **certificates** produced by $\lambda$ and $\lambda'$ are identical.
2. **Automorphism orbits** — the equivalence classes of vertices under the action of the automorphism group.
3. **Automorphism group order** $|\text{Aut}(G)|$.

#### Step 1: Initial Coloring

Group the vertices by their `data` value: vertices with the same `data` belong to the same **cell**. For example, all `Head(42)` vertices are in one cell and all `Slot(7)` vertices in another.

This defines an **ordered partition** $\pi_0 = (C_1, C_2, \dots, C_k)$, where each $C_i$ is a set of vertices.

#### Step 2: 1-WL Color Refinement

The refinement process iteratively splits every non-trivial cell ($|C| > 1$) according to the **neighbour signature** of its vertices with respect to the other cells:

For the signature of a vertex $v \in C_i$ with respect to a target cell $C_j$:

$$
\sigma(v, C_j) = \text{sorted}\bigl[\,(\text{edge\_data}, \text{direction}) \;\big|\; \text{edge } v \to w,\; w \in C_j\,\bigr]
$$

where `direction` is encoded as: `0` = directed outgoing edge, `1` = directed incoming edge, `2` = undirected edge.

If vertices in the same cell have different signatures, the cell is split into several sub-cells by signature. This process repeats until the partition is **equitable** — i.e. for any two vertices $u, v$ in the same cell and any cell $C_j$, $\sigma(u, C_j) = \sigma(v, C_j)$.

**Complexity**: each refinement round produces at most $O(n)$ new cells, and each step takes $O(n^2)$ time to check all (cell, target) pairs. In practice, refinement typically converges within $O(n \log n)$ steps.

#### Step 3: The Individualization–Refinement Search (DFS)

If the partition is still non-discrete after refinement (some cell has $|C| > 1$), the algorithm enters the search phase:

1. **Choose the smallest non-trivial cell** $C_{\min}$.
2. **Individualize** each vertex $v$ of $C_{\min}$: move $v$ out of $C_{\min}$ into a new singleton cell (placed before $C_{\min}$). This breaks the symmetry between $v$ and the other vertices.
3. Re-refine the partition after individualization.
4. Recurse until the partition is discrete (every cell has exactly one vertex), producing a **discrete labeling**.

The search uses an **iterative (stack-based) DFS**; `SearchFrame` records the partition, the path invariant, and whether the node is on the leftmost path at each level.

#### Step 4: Path Invariants and Automorphism Pruning

##### Path Invariants

Every search node records a **path invariant** $\mathcal{I}$ — the sequence of cell sizes of the partitions along the path from the root to the current node:

$$
\mathcal{I} = \bigl[\,|C_1^{(0)}|, \dots, |C_k^{(0)}|\,\bigr] \;\to\; \bigl[\,|C_1^{(1)}|, \dots, |C_{k'}^{(1)}|\,\bigr] \;\to\; \cdots
$$

If two search nodes have different path invariants, they cannot produce isomorphic discrete labelings, so they can be pruned safely. (The current implementation computes and stores this invariant in `SearchFrame`, but invariant-based pruning is not yet active — the only pruning actually performed is the automorphism-orbit pruning.)

##### Automorphism-Orbit Pruning

When the algorithm discovers that two different discrete labelings produce the same certificate, the mapping between them is an automorphism. This automorphism is recorded as an **orbit generator**.

In subsequent search, if the vertex $v$ to be individualized lies in the same orbit as an already processed vertex $w$ under the automorphism group (i.e. there is an automorphism mapping $v$ to $w$), then $v$ is skipped — because individualizing $v$ necessarily produces a subtree isomorphic to the one from individualizing $w$.

Orbits are computed by BFS closure: start from a seed vertex and repeatedly apply the generators (and their inverses) until closure.

#### Step 5: The Canonical Form

During the search, the algorithm maintains the **lexicographically largest certificate** found so far. Whenever a discrete node is reached, its certificate is computed:

$$
\text{cert}(\lambda) = \text{sorted}\bigl[\,(i, j, \text{edge\_data}, \text{dir}) \;\big|\; \text{edge } \lambda^{-1}(i) \to \lambda^{-1}(j)\,\bigr]
$$

If the new certificate is larger, the best labeling is updated; if it is equal, an automorphism generator is recorded.

Final output:

- **Canonical graph**: the graph relabeled by the best labeling.
- **Automorphism group order**: computed by enumerating all permutations generated by the generators via BFS closure (bounded by $n!$, in practice far smaller).
- **Orbit partition**: computed from the generator set by `compute_orbits`.

#### A Complete Example

Consider the expression $T_{\mu\nu} U^{\nu\rho}$, where $T$ and $U$ have no symmetries:

**Graph encoding**:

```
Vertices: H_T (Head), S_mu (Slot), S_nu_T (Slot), H_U (Head), S_nu_U (Slot), S_rho (Slot)
Edges:
  H_T → S_mu    (HeadToSlot(0, 1))
  H_T → S_nu_T  (HeadToSlot(1, 1))
  H_U → S_nu_U  (HeadToSlot(0, 1))
  H_U → S_rho   (HeadToSlot(1, 1))
  S_nu_T — S_nu_U (Contraction(0))   // undirected
```

**Initial coloring**: vertices are grouped by their `data` — `{H_T}` and `{H_U}` form two separate Head cells (different hashes of the function names); slot vertices are grouped by the hash of their labels: `{S_mu}`, `{S_nu_T, S_nu_U}` (both labelled `nu`), `{S_rho}`.

**Refinement**: the non-trivial cell `{S_nu_T, S_nu_U}` has different neighbour signatures with respect to the target cell `{H_T}` ($S_{\nu_T}$ has an incoming edge from $H_T$ — the head-to-slot edge runs from $H_T$ to $S_{\nu_T}$ — while $S_{\nu_U}$ does not), so it is split.

**Result**: after refinement a discrete partition is reached and no search is needed. The canonical form determines the order of the tensor factors and the naming of the dummy indices.

---

### Young Tableaux and Symmetrization Projectors

Young tableaux are the classical tool for representing the irreducible representations of the symmetric group; in tensor algebra they construct tensor components with prescribed symmetries.

#### Definition of a Young Tableau

A **Young tableau** is defined by an integer partition $\lambda = (\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_k > 0)$, where $\lambda_i$ is the number of boxes in the $i$-th row. The total number of boxes $|\lambda| = \sum \lambda_i$ equals the rank of the tensor.

For example, the Young diagram for $\lambda = (2, 1)$:

```
□ □
□
```

#### The Young Symmetrizer

Given a Young tableau $\lambda$, its **Young symmetrizer** $e_\lambda$ acts on a tensor $T$ to produce the component with the symmetry prescribed by $\lambda$. It is defined as follows:

1. **Row symmetrization**: fully symmetrize the indices within each row.
2. **Column antisymmetrization**: fully antisymmetrize the indices within each column.
3. **Composition**: antisymmetrize over columns first, then symmetrize over rows.

Formally, letting $R$ be the direct product of the row symmetry groups and $C$ the direct product of the column antisymmetry groups:

$$
e_\lambda = \left(\sum_{\sigma \in R} \sigma\right) \cdot \left(\sum_{\tau \in C} \operatorname{sgn}(\tau) \cdot \tau\right)
$$

The projection onto a tensor $T$ is:

$$
(e_\lambda \cdot T)_{i_1 \dots i_n} = \sum_{\sigma \in R} \sum_{\tau \in C} \operatorname{sgn}(\tau) \cdot T_{i_{\sigma(\tau(1))} \dots i_{\sigma(\tau(n))}}
$$

#### Computing the Sign of a Permutation

For a given permutation $\pi$, the sign computation of the Young projection proceeds in two steps:

1. **Row-constraint check**: for each position $p$, check whether $\pi(p)$ lies in the same row as $p$. If not (and the row length > 1), then $\operatorname{sgn} = 0$ (the permutation does not preserve the shape of the tableau).

2. **Column parity**: for each column $c$, extract the sub-permutation of $\pi$ restricted to the positions of that column and compute its parity (via cycle decomposition: a cycle of even length contributes a factor $-1$). The parities of all columns multiply to the total sign.

$$
\operatorname{sgn}_\lambda(\pi) = \begin{cases} \displaystyle\prod_{c} (-1)^{\#\text{even cycles in column } c} & \text{if } \pi \text{ preserves rows} \\ 0 & \text{otherwise} \end{cases}
$$

#### The oCAS Implementation: Explicit Expansion

The implementation in `young.rs` uses an **explicit permutation expansion** (not a BSGS group-theoretic method):

1. Expand the row group $R$ (permutations within each row) and the column group $C$ (permutations within each column) into full position maps, for a total of $\prod r_i! \cdot \prod c_j!$ combinations.
2. For each combination $(\sigma, \tau)$, form the composed permutation $\tau \circ \sigma$; the term coefficient is the column permutation parity $\operatorname{sgn}(\tau)$.
3. Each term contributes $\operatorname{sgn}(\tau) \cdot T(i_{\tau\sigma(1)}, \dots, i_{\tau\sigma(n)})$.
4. Sum all terms — this is the classical Young symmetrizer $c_\lambda = a_\lambda b_\lambda$, with **no** normalization by the product of hook lengths; every term has coefficient $\pm 1$ only. The result is therefore a projector up to the scalar $\kappa = (\prod r_i! \cdot \prod c_j!) / \dim(\lambda)$ ($c_\lambda^2 = \kappa \cdot c_\lambda$).

**Fully antisymmetric tableau** $\lambda = (1, 1, \dots, 1)$: the row group is trivial and the column group is the whole $S_n$, producing the standard alternating sum.

**Fully symmetric tableau** $\lambda = (n)$: one element per column, so column antisymmetrization degenerates to the identity and row symmetrization produces the sum over all permutations.

---

### Symmetry Specifications

The `SymmetrySpec` struct declares the slot symmetries of a tensor function head; it is the bridge between graph encoding and Young projection.

#### Three Types of Symmetry

| Type | Field | Graph-encoding behavior | Young-projection behavior |
|---|---|---|---|
| Symmetric subsets | `symmetric_subsets` | slots are marked hidden (`is_slot_hidden = true`); the graph isomorphism may permute them freely | preserves the corresponding row symmetry |
| Antisymmetric subsets | `antisymmetric_subsets` | slots are visible and encoded exactly like ordinary slots (no special marker; the antisymmetry signs are handled at the expression level by the Young projection) | preserves the corresponding column antisymmetry |
| Cyclic permutations | `cyclic` | slots are marked hidden | symmetrization over cyclic permutations |

#### The `is_slot_hidden` Test

A slot $p$ is marked hidden if and only if:

$$
\text{is\_hidden}(p) \iff \exists S \in \text{symmetric\_subsets}: p \in S \;\lor\; \exists C = \text{cyclic}: p \in C
$$

Hidden slots use the uniform color `Slot(0)` in the graph, letting the graph-isomorphism engine permute them freely, which guarantees $T(a, b) = T(b, a)$ for symmetric tensors.

#### TensorRegistry

`TensorRegistry` manages the specifications and index groups of all tensors:

- `register(name, spec)`: register a tensor function head together with its symmetry specification.
- `set_index_group(label, group)`: set the group to which an index label belongs (e.g. spacetime indices group 1, internal indices group 2).
- `index_group(label)`: query the index group (0 = default/unassigned).

The role of index groups is to prevent cross-dimension conflicts when renaming dummy indices. For example, a spacetime dummy index $\mu$ and an internal dummy index $i$ use different namespaces `d1_0`, `d1_1`, … and `d2_0`, `d2_1`, ….

---

## Implementation in oCAS

### Module Structure

```
ocas-atom/src/tensor/
├── canon.rs      # graph encoding + canonicalization entry
├── graph.rs      # McKay refinement-individualization engine (about 1100 lines)
├── young.rs      # Young tableau projection
├── spec.rs       # SymmetrySpec + TensorRegistry
└── dummy.rs      # dummy-index management
```

### `canon.rs`: Graph Encoding and Canonicalization

**Entry function** `canonicalize_tensors(ctx, expr, registry)`:

1. If the expression is a sum (`Add`), canonicalize each term and verify that the free indices are consistent.
2. A monomial goes to `canonicalize_single_term`:
   - **Fast path**: if the tensor has symmetric subsets, no antisymmetric subsets, and all slots are hidden (`is_slot_hidden` is true for all), simply sort the arguments by the alphabetical order of their labels; no graph isomorphism is needed.
   - **General path**: encode the graph with `tensor_to_graph`, compute the canonical form with `g.canonize()`, and rebuild the expression with `reconstruct`.

**Encoding logic of `tensor_to_graph`**:

- For each factor `Fun(name, args)` in the product:
  - create a head vertex `Head(hash(name))`.
  - pre-sort the symmetric slots by the alphabetical order of their labels, so that $T(a,b)$ and $T(b,a)$ generate the same graph.
  - create a slot vertex for each argument; hidden slots use `Slot(0)` and visible slots use `Slot(hash(label))`.
  - add the directed edge `Head → Slot`.
- For pairs of slots whose label occurs exactly twice in the product, add the undirected contraction edge.

**Rebuilding logic of `reconstruct`**:

- Walk the head vertices of the canonical graph and collect the tensor factors in canonical order.
- For each head's outgoing edges (sorted by `orig_pos` or label order), determine the slot labels:
  - contracted slots: assign the canonical dummy names `d0`, `d1`, …
  - free slots: keep the original labels.
- Assemble the factors into a product or a sum expression.

### `graph.rs`: The McKay Algorithm

This is the largest and most self-contained algorithmic module in oCAS (about 1100 lines), implementing a complete nauty-style graph-isomorphism engine.

**Data structures**:

- `Graph<N, H, E>`: a three-parameter generic graph. `N` = visible vertex data (participates in comparison), `H` = hidden vertex data (does not participate), `E` = edge data (participates).
- `Partition`: an ordered partition; `cells[i]` is the vertex list of the $i$-th cell and `cell_of[v]` is the index of the cell containing vertex $v$.
- `SearchFrame`: a DFS search stack frame holding the partition, the path invariant, and whether the node is on the leftmost path.
- `CanonicalForm<N, H, E>`: the output result, containing `vertex_map` (original → canonical mapping), `orbits` (automorphism orbits), `automorphism_group_size` (group order), and `graph` (the canonical graph).

**Key implementation details**:

- **Refinement** (`Partition::refine`): uses the `stable_below` pointer optimization — cells that are already stable no longer participate in split checks. When a cell is split, `stable_below` retreats to that cell's index, ensuring subsequent refinement covers all affected cells.
- **Neighbour signatures** (`cell_signatures` / `vertex_signature`): for each vertex $v$, collect the `(edge_data, direction)` pairs to all neighbours in the target cell $C_j$, and sort them into a signature. `direction` distinguishes directed outgoing (0), directed incoming (1), and undirected (2) edges.
- **Certificate** (`Graph::certificate`): convert all edges of the graph into a list of `(i, j, data, dir)` quadruples under the canonical labeling, sorted to serve as the key for lexicographic comparison.
- **Automorphism group order** (`group_size`): enumerates all permutations generated by the generators via BFS closure (**full enumeration**, no truncation; bounded by $n!$, in practice far smaller).
- **Orbit computation** (`compute_orbits`): BFS starts from each unvisited vertex and expands along the generators and their inverses until closure.

### `young.rs`: Young Projection

**The `YoungTableau` struct**:

- `row_lengths: Vec<usize>`: the number of boxes in each row; e.g. `[2, 1]` denotes $\lambda = (2,1)$.
- `total_boxes()`: the total number of boxes.
- `sign_of_permutation(perm)`: compute the sign of a permutation under the given tableau ($-1, 0, +1$).

**`young_project(ctx, tensor_expr, tableau)`**:

- only handles expressions of the form `Fun(name, args)`.
- generates all permutations with Heap's algorithm ($O(n!)$ time, constant space).
- computes the sign of each permutation and, when non-zero, constructs the signed term.
- returns the sum of all terms (or a single term, or zero).

**Note**: this is an explicit expansion, suitable for low-rank tensors ($n \lesssim 8$). For high-rank tensors, a BSGS (basic Schreier–Sims) method should be considered, but it is not currently implemented in oCAS.

**Known limitation**: the column-parity code of `sign_of_permutation` panics (via `unwrap`) on permutations that **preserve rows but move elements across columns** (the column-position lookup for the element's new position fails). The fully symmetric tableau $(n)$ and the fully antisymmetric tableau $(1, \ldots, 1)$ never trigger this; general shapes (such as $(2,1)$) can panic during the explicit expansion, so care is needed.

### `spec.rs`: Symmetry Specifications

**The three constructors of `SymmetrySpec`**:

- `none()`: no symmetry.
- `fully_symmetric(rank)`: all slots fully symmetric.
- `fully_antisymmetric(rank)`: all slots fully antisymmetric.

**`is_slot_hidden(pos)`**: checks whether a slot belongs to a symmetric subset or a cyclic subset. Hidden slots do not participate in the graph-canonicalization comparison, allowing the graph-isomorphism engine to permute them freely.

### `dummy.rs`: Dummy-Index Management

**`refresh_dummies(ctx, expr, registry)`**:

1. Collect the occurrence counts of all `Fun` node arguments in the expression.
2. Identify the dummy indices (occurring exactly twice).
3. For each dummy index, query its index group (`registry.index_group`) and assign a number within the group.
4. Rename: dummy indices of group 0 are named `d0`, `d1`, …; those of group $g$ are named `d{g}_0`, `d{g}_1`, ….

**Error variants** (declared in `dummy.rs`, but not raised by `refresh_dummies`):

- `OverContracted(Symbol)`: an index occurs more than twice.
- `BadContraction(Symbol)`: the two occurrences have the same variance (not an upper/lower pair).

(In the canonicalization path, a label occurring more than twice raises `TensorCanonError::ContractedMoreThanOnce` from `canon.rs`.)

---

## Advanced Topics

### Hidden Slots and the Fast Path

For fully symmetric tensors (such as the metric $g_{\mu\nu}$), all slots are marked hidden. All slot vertices in the graph then share the same color `Slot(0)`, and the automorphism group of the McKay algorithm becomes very large ($|S_n| = n!$).

To avoid this cost, `canonicalize_single_term` takes the **fast path** when it detects a fully symmetric tensor: sort the arguments directly by the alphabetical order of their labels without building a graph. This reduces the canonicalization complexity from $O(n!)$ to $O(n \log n)$.

### Uses of Automorphism-Group Information

The McKay algorithm outputs not only the canonical form, but also the orbits and the order of the automorphism group. This information has important applications in tensor algebra:

- **Orbits**: identify which slots are equivalent under the symmetries, useful for simplifying sums.
- **Group order**: $|\text{Aut}(G)|$ is a measure of the intrinsic symmetry of the expression, useful for estimating the size of a symmetrization expansion.

### Relation to nauty

`graph.rs` is a standalone nauty-style implementation with no external dependencies. It implements the core components of McKay's algorithm (1981/2014):

- 1-WL color refinement (equivalent to nauty's `refine`)
- individualization–refinement DFS (equivalent to `search`)
- path-invariant recording (equivalent to nauty's `compare`; pruning based on it is not yet active)
- automorphism-orbit pruning (equivalent to nauty's `orbits`)

The main differences from the original nauty:

- uses adjacency lists instead of bit matrices, suitable for sparse graphs.
- the certificate uses a sorted edge list instead of a canonical adjacency matrix.
- the hidden data (the `H` parameter) allows storing metadata that does not participate in comparison.

### Extending to Higher Symmetries

For more complex symmetries (such as the algebraic Bianchi identities of the Riemann curvature tensor $R_{abcd}$), the current system approximates them by combining `symmetric_subsets` and `antisymmetric_subsets`. Full Young symmetrization is provided by the `young_project` function, but it produces an expanded sum rather than a compact representation.

---

## References

1. **Cvitanović, P.** *Group Theory: Birdtracks, Lie's, and Exceptional Groups.* Princeton University Press, 2008. — The classic reference on tensor invariants and representations of symmetry groups.
2. **McKay, B. D.** *Practical Graph Isomorphism.* Congressus Numerantium, 30:45–87, 1981 (updated 2014). — The original paper on the nauty algorithm.
3. **Faugère, J.-C.** *A new efficient algorithm for computing Gröbner bases (F4).* Journal of Pure and Applied Algebra, 139(1–3):61–88, 1999. — The F4 algorithm, related to the matrix method for tensor contractions.
4. **Cox, D., Little, J., O'Shea, D.** *Ideals, Varieties, and Algorithms.* Springer, 4th ed., 2015. — Polynomial-algebra foundations; reference for Gröbner basis theory.
5. **Hamermesh, M.** *Group Theory and Its Application to Physical Problems.* Dover, 1989. — Physical applications of Young tableaux and tensor symmetries.
