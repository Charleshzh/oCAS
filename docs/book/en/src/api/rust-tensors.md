# Rust API Reference: Tensors

oCAS provides basic tensor algebra with explicit index management. A `Tensor` carries named index slots (upper/lower), optional symmetry metadata, and can be converted into an `Atom` expression for symbolic processing. Since 0.22.0, a graph-isomorphism-based canonicalization engine and Young projectors are also available.

**Module paths**: `ocas_atom::tensor` (`mod.rs`), `ocas_atom::tensor::spec`, `ocas_atom::tensor::canon`, `ocas_atom::tensor::dummy`, `ocas_atom::tensor::young`

**Import**:

```rust
use ocas_atom::tensor::{
    IndexPosition, IndexSlot, Symmetry, Tensor,
    contract, Contracted, TensorProduct,
    symmetrise_sign,
};
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::{canonicalize_tensors, CanonicalTensor, TensorCanonError};
use ocas_atom::tensor::dummy::{refresh_dummies, DummyError};
use ocas_atom::tensor::young::{YoungTableau, young_project};
```

---

## IndexPosition

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IndexPosition {
    Upper,
    Lower,
}
```

**Description**: The position (variance) of an index on a tensor. `Upper` denotes a contravariant index (superscript), `Lower` a covariant index (subscript). Contraction requires two indices with the same label to have opposite positions.

**Variants**:

| Variant | Meaning |
|---|---|
| `Upper` | Contravariant index (superscript); corresponds to a contravariant vector in physics |
| `Lower` | Covariant index (subscript); corresponds to a covariant vector in physics |

**Example**:

```rust
use ocas_atom::tensor::IndexPosition;

assert_eq!(IndexPosition::Upper, IndexPosition::Upper);
assert_ne!(IndexPosition::Upper, IndexPosition::Lower);
// Output: all assertions pass
```

---

## IndexSlot

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IndexSlot<'a> { /* private fields */ }
```

**Description**: A single index slot of a tensor, consisting of an index label (`Atom`) and a position (`IndexPosition`). It is a `Copy` type and can be copied cheaply.

### IndexSlot::new

**Signature**: `pub fn new(label: Atom<'a>, position: IndexPosition) -> Self`

**Description**: Creates a new index slot.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `label` | `Atom<'a>` | Index label expression (usually a variable name such as `mu` or `i`) |
| `position` | `IndexPosition` | Index position: `Upper` or `Lower` |

**Returns**: `IndexSlot<'a>`

**Example**:

```rust
use ocas_atom::tensor::{IndexSlot, IndexPosition};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mu = IndexSlot::new(ctx.var("mu"), IndexPosition::Upper);
let nu = IndexSlot::new(ctx.var("nu"), IndexPosition::Lower);
assert_eq!(mu.position(), IndexPosition::Upper);
assert_eq!(nu.position(), IndexPosition::Lower);
// Output: all assertions pass
```

### IndexSlot::label

**Signature**: `pub fn label(&self) -> Atom<'a>`

**Description**: Returns the index label expression.

**Returns**: `Atom<'a>` — the `Atom` handle of the index label.

### IndexSlot::position

**Signature**: `pub fn position(&self) -> IndexPosition`

**Description**: Returns the position of the index (`Upper` or `Lower`).

**Returns**: `IndexPosition`

**See also**: [IndexPosition](#indexposition)

---

## Symmetry

**Signature**:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Symmetry {
    None,
    Symmetric,
    Antisymmetric,
}
```

**Description**: The symmetry of a tensor's index slots. This is **advisory** metadata — operations such as `contract` do not automatically symmetrize; the symmetry is used by downstream consumers (such as the canonicalization engine or Young projectors).

**Variants**:

| Variant | Meaning |
|---|---|
| `None` | No symmetry (generic tensor) |
| `Symmetric` | Symmetric: swapping any two index slots leaves the tensor unchanged |
| `Antisymmetric` | Antisymmetric: swapping any two index slots flips the sign |

**Example**:

```rust
use ocas_atom::tensor::Symmetry;

assert_eq!(Symmetry::None, Symmetry::None);
assert_ne!(Symmetry::Symmetric, Symmetry::Antisymmetric);
// Output: all assertions pass
```

---

## Tensor

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tensor<'a> { /* private fields */ }
```

**Description**: A named tensor object, consisting of a name (`Symbol`), a list of index slots (`Vec<IndexSlot>`), and a symmetry (`Symmetry`). It can be lowered to a standard `Atom` expression node via `to_atom`.

### Tensor::new

**Signature**: `pub fn new(name: Symbol, slots: Vec<IndexSlot<'a>>) -> Self`

**Description**: Creates a new tensor; the default symmetry is `Symmetry::None`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `Symbol` | Tensor name (interned string, e.g. `"g"`, `"Riemann"`) |
| `slots` | `Vec<IndexSlot<'a>>` | List of index slots; the order is meaningful |

**Returns**: `Tensor<'a>`

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let slots = vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
];
let t = Tensor::new(Symbol::new("g"), slots);
assert_eq!(t.rank(), 2);
assert_eq!(t.name().as_str(), "g");
// Output: all assertions pass
```

### Tensor::with_symmetry

**Signature**: `pub fn with_symmetry(mut self, symmetry: Symmetry) -> Self`

**Description**: Builder method that sets the tensor's symmetry. Consumes and returns `self`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `symmetry` | `Symmetry` | The desired symmetry |

**Returns**: `Self` — the tensor with the symmetry set.

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, Symmetry};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let slots = vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Upper),
];
let sym = Tensor::new(Symbol::new("g"), slots.clone())
    .with_symmetry(Symmetry::Symmetric);
let anti = Tensor::new(Symbol::new("epsilon"), slots)
    .with_symmetry(Symmetry::Antisymmetric);

assert_eq!(sym.symmetry(), Symmetry::Symmetric);
assert_eq!(anti.symmetry(), Symmetry::Antisymmetric);
// Output: all assertions pass
```

### Tensor::name

**Signature**: `pub fn name(&self) -> Symbol`

**Description**: Returns the tensor name.

**Returns**: `Symbol` — the interned string handle.

### Tensor::slots

**Signature**: `pub fn slots(&self) -> &[IndexSlot<'a>]`

**Description**: Returns a slice reference to the index slots.

**Returns**: `&[IndexSlot<'a>]`

### Tensor::symmetry

**Signature**: `pub fn symmetry(&self) -> Symmetry`

**Description**: Returns the tensor's symmetry.

**Returns**: `Symmetry`

### Tensor::rank

**Signature**: `pub fn rank(&self) -> usize`

**Description**: Returns the rank of the tensor (the number of index slots).

**Returns**: `usize`

### Tensor::dummy_labels

**Signature**: `pub fn dummy_labels(&self) -> Vec<Atom<'a>>`

**Description**: Returns the list of dummy index labels. A dummy index is a label that occurs exactly twice across all slots (once `Upper` and once `Lower`), i.e. a contraction candidate.

**Returns**: `Vec<Atom<'a>>` — the list of dummy index labels.

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// A^μ_μ: mu occurs twice (Upper + Lower) → dummy index
let t = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("mu"), IndexPosition::Lower),
]);
assert_eq!(t.dummy_labels().len(), 1);
// Output: 1 dummy index
```

### Tensor::to_atom

**Signature**: `pub fn to_atom(&self, ctx: &'a AtomArena<'a>) -> Atom<'a>`

**Description**: Lowers the tensor to a standard `Atom` function node `name(slot₁, slot₂, …)`. No symmetrization is applied — the atom preserves `self`'s slot order.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |

**Returns**: `Atom<'a>` — the function node expression.

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let t = Tensor::new(Symbol::new("g"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
]);
let atom = t.to_atom(&ctx);
println!("{}", atom);  // g(mu, nu)
// Output: g(mu, nu)
```

**See also**: [Atom](./rust-expressions.md#atom)

---

## contract

**Signature**:

```rust
pub fn contract<'a>(
    ctx: &'a AtomArena<'a>,
    a: &Tensor<'a>,
    b: &Tensor<'a>,
) -> Contracted<'a>
```

**Description**: Contracts two tensors by summing over shared dummy indices. Two slots with the same label but opposite positions are contracted. The result retains the surviving free indices, in the order `(a's free slots, b's free slots)`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `a` | `&Tensor<'a>` | The first tensor |
| `b` | `&Tensor<'a>` | The second tensor |

**Returns**: `Contracted<'a>`
- With no contraction pair, returns `Contracted::Product` containing the two original tensors
- With a partial contraction, returns `Contracted::Product` containing a new tensor carrying the free slots
- With a full contraction (no free slots), returns `Contracted::Scalar`

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, contract, Contracted};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// A^μ B_μ → scalar
let a = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
]);
let b = Tensor::new(Symbol::new("B"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Lower),
]);
match contract(&ctx, &a, &b) {
    Contracted::Scalar(expr) => println!("scalar: {}", expr),
    Contracted::Product(tp) => {
        for f in &tp.factors {
            println!("free: {} (rank {})", f.name().as_str(), f.rank());
        }
    }
}
// Output: scalar: A(mu)*B(mu)
```

**See also**: [Contracted](#contracted), [TensorProduct](#tensorproduct)

---

## Contracted

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Contracted<'a> {
    Product(TensorProduct<'a>),
    Scalar(Atom<'a>),
}
```

**Description**: The result of contracting two tensors. When no free indices survive it is a scalar expression; otherwise it is a tensor product with free slots.

**Variants**:

| Variant | Type | Description |
|---|---|---|
| `Product` | `TensorProduct<'a>` | Partial or no contraction; a tensor product with free indices |
| `Scalar` | `Atom<'a>` | Full contraction to a scalar expression |

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, contract, Contracted};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// A^μ_ν B^ν_ρ → partial contraction, leaving A^μ_ρ
let a = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
]);
let b = Tensor::new(Symbol::new("B"), vec![
    IndexSlot::new(ctx.var("nu"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("rho"), IndexPosition::Lower),
]);
match contract(&ctx, &a, &b) {
    Contracted::Product(tp) => {
        assert_eq!(tp.factors.len(), 1);
        assert_eq!(tp.factors[0].rank(), 2); // mu, rho
    }
    Contracted::Scalar(_) => panic!("expected partial contraction"),
}
// Output: partial contraction, 2 free indices remain
```

---

## TensorProduct

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TensorProduct<'a> {
    pub factors: Vec<Tensor<'a>>,
}
```

**Description**: A tensor product retaining the free slots after contraction. `factors` contains the tensors that survive the contraction, with slots concatenated in the order `(a's free slots, b's free slots)`.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `factors` | `Vec<Tensor<'a>>` | The list of tensor factors after contraction |

**See also**: [contract](#contract), [Contracted](#contracted)

---

## symmetrise_sign

**Signature**: `pub fn symmetrise_sign(tensor: &Tensor<'_>) -> i64`

**Description**: Sorts the tensor's index slots by label in ascending order and returns the sign of the antisymmetrization.

- `Symmetry::None` or `Symmetry::Symmetric`: returns `+1`
- `Symmetry::Antisymmetric`: returns the permutation parity (even permutation `+1`, odd permutation `-1`)

**Note**: This is not full canonicalization under the permutation group (that requires graph isomorphism); it merely provides a stable order for comparing equivalent symmetric tensors.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `tensor` | `&Tensor<'_>` | The input tensor |

**Returns**: `i64` — `+1` or `-1`.

**Example**:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, Symmetry, symmetrise_sign};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let slots = vec![
    IndexSlot::new(ctx.var("a"), IndexPosition::Upper),
    IndexSlot::new(ctx.var("b"), IndexPosition::Upper),
];
let anti = Tensor::new(Symbol::new("epsilon"), slots)
    .with_symmetry(Symmetry::Antisymmetric);

let sign = symmetrise_sign(&anti);
// If a < b the permutation is the identity, sign = +1
// Output: +1
```

---

## SymmetrySpec

**Signature**:

```rust
#[derive(Debug, Clone, Default)]
pub struct SymmetrySpec {
    pub symmetric_subsets: Vec<Vec<usize>>,
    pub antisymmetric_subsets: Vec<Vec<usize>>,
    pub cyclic: Option<Vec<usize>>,
}
```

**Description**: A fine-grained symmetry specification for tensors, used by the canonicalization engine. It is more flexible than the `Symmetry` enum — different subsets of slot positions can be assigned different symmetry behavior.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `symmetric_subsets` | `Vec<Vec<usize>>` | List of symmetric slot subsets. Slots within the same subset are interchangeable |
| `antisymmetric_subsets` | `Vec<Vec<usize>>` | List of antisymmetric slot subsets. Swapping two slots within a subset flips the sign |
| `cyclic` | `Option<Vec<usize>>` | A slot subset undergoing cyclic permutations |

### SymmetrySpec::none

**Signature**: `pub fn none() -> Self`

**Description**: No symmetry — all slots are independent.

**Returns**: `SymmetrySpec`

**Example**:

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::none();
assert!(spec.symmetric_subsets.is_empty());
assert!(spec.antisymmetric_subsets.is_empty());
assert!(spec.cyclic.is_none());
// Output: all assertions pass
```

### SymmetrySpec::fully_symmetric

**Signature**: `pub fn fully_symmetric(rank: usize) -> Self`

**Description**: All slots are fully symmetric. Creates a symmetric subset containing all slots `0..rank`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rank` | `usize` | The rank of the tensor (total number of slots) |

**Returns**: `SymmetrySpec`

**Example**:

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::fully_symmetric(3);
assert_eq!(spec.symmetric_subsets, vec![vec![0, 1, 2]]);
// Output: symmetric subset = [[0, 1, 2]]
```

### SymmetrySpec::fully_antisymmetric

**Signature**: `pub fn fully_antisymmetric(rank: usize) -> Self`

**Description**: All slots are fully antisymmetric. Creates an antisymmetric subset containing all slots `0..rank`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `rank` | `usize` | The rank of the tensor (total number of slots) |

**Returns**: `SymmetrySpec`

**Example**:

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::fully_antisymmetric(4);
assert_eq!(spec.antisymmetric_subsets, vec![vec![0, 1, 2, 3]]);
// Output: antisymmetric subset = [[0, 1, 2, 3]]
```

### SymmetrySpec::is_slot_hidden

**Signature**: `pub fn is_slot_hidden(&self, pos: usize) -> bool`

**Description**: Checks whether the given slot position should be "hidden" in the graph encoding (i.e. it does not participate in canonicalization comparisons). Symmetric and cyclic slots are hidden.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `pos` | `usize` | The slot position index |

**Returns**: `bool` — `true` if the slot is hidden.

**Example**:

```rust
use ocas_atom::tensor::spec::SymmetrySpec;

let spec = SymmetrySpec::fully_symmetric(2);
assert!(spec.is_slot_hidden(0));   // in symmetric subset
assert!(spec.is_slot_hidden(1));   // in symmetric subset

let none = SymmetrySpec::none();
assert!(!none.is_slot_hidden(0));  // not in any subset
// Output: all assertions pass
```

---

## TensorRegistry

**Signature**:

```rust
#[derive(Debug, Clone, Default)]
pub struct TensorRegistry { /* private fields */ }
```

**Description**: The complete registry for tensor canonicalization. It records which function heads are tensors and their symmetry specifications, and also manages the assignment of index dimension groups (preventing dummy-index renaming conflicts across dimensions).

### TensorRegistry::new

**Signature**: `pub fn new() -> Self`

**Description**: Creates an empty registry.

**Returns**: `TensorRegistry`

### TensorRegistry::register

**Signature**: `pub fn register(&mut self, name: Symbol, spec: SymmetrySpec)`

**Description**: Registers a tensor name together with its symmetry specification.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `Symbol` | The tensor function head name |
| `spec` | `SymmetrySpec` | The symmetry specification |

**Example**:

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::Symbol;

let mut reg = TensorRegistry::new();
reg.register(Symbol::new("g"), SymmetrySpec::fully_symmetric(2));
reg.register(Symbol::new("Riemann"), SymmetrySpec::none());
reg.register(Symbol::new("epsilon"), SymmetrySpec::fully_antisymmetric(4));
// Output: 3 tensors registered
```

### TensorRegistry::set_index_group

**Signature**: `pub fn set_index_group(&mut self, label: Symbol, group: u64)`

**Description**: Sets the dimension-group identifier of an index label. Dummy indices from different groups are never renamed to the same canonical name, avoiding cross-dimension conflicts.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `label` | `Symbol` | The index label (e.g. `"mu"`, `"i"`) |
| `group` | `u64` | The group identifier. `0` means default/ungrouped |

**Example**:

```rust
use ocas_atom::tensor::spec::TensorRegistry;
use ocas_atom::Symbol;

let mut reg = TensorRegistry::new();
reg.set_index_group(Symbol::new("mu"), 1);  // spacetime index
reg.set_index_group(Symbol::new("i"), 2);   // internal index
assert_eq!(reg.index_group(Symbol::new("mu")), 1);
assert_eq!(reg.index_group(Symbol::new("i")), 2);
assert_eq!(reg.index_group(Symbol::new("other")), 0); // default group
// Output: all assertions pass
```

### TensorRegistry::spec

**Signature**: `pub fn spec(&self, name: Symbol) -> Option<&SymmetrySpec>`

**Description**: Looks up the symmetry specification of a registered tensor.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `name` | `Symbol` | The tensor name |

**Returns**: `Option<&SymmetrySpec>` — `None` if not registered.

### TensorRegistry::index_group

**Signature**: `pub fn index_group(&self, label: Symbol) -> u64`

**Description**: Looks up the dimension group of an index label. Labels that have not been set return `0`.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `label` | `Symbol` | The index label |

**Returns**: `u64` — the group identifier.

**See also**: [SymmetrySpec](#symmetryspec)

---

## canonicalize_tensors

**Signature**:

```rust
pub fn canonicalize_tensors<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<CanonicalTensor<'a>, TensorCanonError>
```

**Description**: Canonicalizes a tensor expression to a unique canonical form independent of index naming. Internally, the expression is encoded as a graph (tensor head → head vertex, index slot → slot vertex, contraction → edge), the McKay refinement-individualization algorithm computes canonical labels, and the expression is then rebuilt with renamed dummy indices.

For sum expressions (`Add`), each term is canonicalized individually and all terms are verified to have the same set of free indices.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `expr` | `Atom<'a>` | The tensor expression to canonicalize (`Mul` of `Fun` nodes) |
| `registry` | `&TensorRegistry` | The tensor registry (symmetry specifications + index groups) |

**Returns**: `Result<CanonicalTensor<'a>, TensorCanonError>`

**Errors**:

| Error variant | Description |
|---|---|
| `ContractedMoreThanOnce(Symbol)` | An index label is contracted more than twice |
| `BadContraction(Symbol)` | Two slots with the same label have the same position (not an upper/lower pair) |
| `NotATensor(Symbol)` | A function head not registered in the registry appears in the expression |
| `InconsistentOpenIndices` | The free indices of the terms in a sum are inconsistent |
| `UnsupportedPower` | An unsupported power node appears in the expression |

**Example**:

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mut reg = TensorRegistry::new();
reg.register(Symbol::new("T"), SymmetrySpec::none());

// T(i, j) and T(j, i) are different (no symmetry)
let t_ij = ctx.fun("T", &[ctx.var("i"), ctx.var("j")]);
let t_ji = ctx.fun("T", &[ctx.var("j"), ctx.var("i")]);

let ct1 = canonicalize_tensors(&ctx, t_ij, &reg).unwrap();
let ct2 = canonicalize_tensors(&ctx, t_ji, &reg).unwrap();
// Without symmetry, T(i,j) ≠ T(j,i)
// Note: without symmetry, T(i,j) and T(j,i) have different canonical forms
```

**See also**: [CanonicalTensor](#canonicaltensor), [TensorCanonError](#tensorcanonerror), [TensorRegistry](#tensorregistry)

---

## CanonicalTensor

**Signature**:

```rust
#[derive(Debug, Clone)]
pub struct CanonicalTensor<'a> {
    pub canonical_form: Atom<'a>,
    pub external_indices: Vec<Atom<'a>>,
    pub dummy_indices: Vec<Atom<'a>>,
}
```

**Description**: The result of tensor canonicalization.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `canonical_form` | `Atom<'a>` | The canonicalized expression (dummy indices renamed) |
| `external_indices` | `Vec<Atom<'a>>` | The list of free (external) indices |
| `dummy_indices` | `Vec<Atom<'a>>` | The list of dummy indices (renamed to `d0`, `d1`, …) |

**Example**:

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mut reg = TensorRegistry::new();
reg.register(Symbol::new("T"), SymmetrySpec::none());
reg.register(Symbol::new("U"), SymmetrySpec::none());

let prod = ctx.mul(&[ctx.fun("T", &[ctx.var("i"), ctx.var("j")]),
                      ctx.fun("U", &[ctx.var("j"), ctx.var("k")])]);
let ct = canonicalize_tensors(&ctx, prod, &reg).unwrap();

println!("Canonical form: {}", ct.canonical_form);
println!("External indices: {:?}", ct.external_indices);
println!("Dummy indices: {:?}", ct.dummy_indices);
// Output: canonical form contains renamed dummy indices; external indices are i, k
```

---

## TensorCanonError

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TensorCanonError {
    ContractedMoreThanOnce(Symbol),
    BadContraction(Symbol),
    NotATensor(Symbol),
    InconsistentOpenIndices,
    UnsupportedPower,
}
```

**Description**: Errors that may arise during tensor canonicalization.

**Variants**:

| Variant | Description |
|---|---|
| `ContractedMoreThanOnce(Symbol)` | An index label is contracted more than twice (occurs 3 or more times) |
| `BadContraction(Symbol)` | Two slots with the same label have the same variance (not a valid upper/lower pair) |
| `NotATensor(Symbol)` | A function head not registered in the `TensorRegistry` appears |
| `InconsistentOpenIndices` | The sets of free indices of the terms in a sum are inconsistent |
| `UnsupportedPower` | An unsupported power node appears in the expression |

---

## refresh_dummies

**Signature**:

```rust
pub fn refresh_dummies<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<Atom<'a>, DummyError>
```

**Description**: Renames the dummy indices in a tensor expression to avoid conflicts with external (free) indices. Labels that occur exactly twice (once as a superscript, once as a subscript) are replaced with new names assigned per dimension group (`d0`, `d1`, …). External indices remain unchanged.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `expr` | `Atom<'a>` | The tensor expression to process |
| `registry` | `&TensorRegistry` | The tensor registry (used for dimension-group assignment) |

**Returns**: `Result<Atom<'a>, DummyError>`

**Errors**:

| Error variant | Description |
|---|---|
| `OverContracted(Symbol)` | An index label occurs more than twice |
| `BadContraction(Symbol)` | Two slots with the same label have the same variance |

**Example**:

```rust
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::dummy::refresh_dummies;
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let mut reg = TensorRegistry::new();
reg.register(Symbol::new("T"), SymmetrySpec::none());
reg.register(Symbol::new("U"), SymmetrySpec::none());

// T(i,j) * U(j,i): i and j both occur exactly twice → dummy indices
let expr = ctx.mul(&[ctx.fun("T", &[ctx.var("i"), ctx.var("j")]),
                      ctx.fun("U", &[ctx.var("j"), ctx.var("i")])]);
let refreshed = refresh_dummies(&ctx, expr, &reg).unwrap();
println!("{}", refreshed);
// Output: dummy indices renamed to d0, d1 (assigned per dimension group)
```

**See also**: [DummyError](#dummyerror), [TensorRegistry](#tensorregistry)

---

## DummyError

**Signature**:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DummyError {
    OverContracted(Symbol),
    BadContraction(Symbol),
}
```

**Description**: The error type for dummy-index operations.

**Variants**:

| Variant | Description |
|---|---|
| `OverContracted(Symbol)` | An index label occurs more than twice (not a valid contraction candidate) |
| `BadContraction(Symbol)` | Two slots with the same label have the same variance (not an upper/lower pair) |

---

## YoungTableau

**Signature**:

```rust
#[derive(Debug, Clone)]
pub struct YoungTableau {
    pub row_lengths: Vec<usize>,
}
```

**Description**: A Young tableau, whose shape is defined by a list of row lengths (e.g. `[2, 1]` denotes □□/□). The Young projector implements the classical symmetrizer $c_\lambda = a_\lambda \cdot b_\lambda$: it sums over all combinations of row permutations $\sigma \in R$ and column permutations $\tau \in C$, each term carrying the column permutation parity $\operatorname{sgn}(\tau)$.

This is an **explicit** expansion (not a BSGS group-theoretic implementation): the result is a sum of $\prod r_i! \cdot \prod c_j!$ terms (products of row/column factorials), each with sign ±1 and no hook-length normalization.

**Fields**:

| Field | Type | Description |
|---|---|---|
| `row_lengths` | `Vec<usize>` | The number of boxes in each row |

### YoungTableau::new

**Signature**: `pub fn new(row_lengths: Vec<usize>) -> Self`

**Description**: Creates a Young tableau from row lengths.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `row_lengths` | `Vec<usize>` | The list of boxes per row. The total number of boxes should equal the rank of the tensor |

**Returns**: `YoungTableau`

**Example**:

```rust
use ocas_atom::tensor::young::YoungTableau;

// [2, 1] = □□/□ (3 boxes)
let tab = YoungTableau::new(vec![2, 1]);
assert_eq!(tab.total_boxes(), 3);

// [1, 1, 1] = fully antisymmetric
let anti = YoungTableau::new(vec![1, 1, 1]);
assert_eq!(anti.total_boxes(), 3);

// [3] = fully symmetric
let sym = YoungTableau::new(vec![3]);
assert_eq!(sym.total_boxes(), 3);
// Output: all assertions pass
```

### YoungTableau::total_boxes

**Signature**: `pub fn total_boxes(&self) -> usize`

**Description**: Returns the total number of boxes in the tableau (equal to the rank of the tensor).

**Returns**: `usize`

---

## young_project

**Signature**:

```rust
pub fn young_project<'a>(
    ctx: &'a AtomArena<'a>,
    tensor_expr: Atom<'a>,
    tableau: &YoungTableau,
) -> Atom<'a>
```

**Description**: Applies the Young projector to a tensor expression. Expands $T(i_1, i_2, \dots, i_n)$ into a sum over permutations:

$$\sum_\sigma \text{sign}(\sigma) \cdot T(i_{\sigma(1)}, i_{\sigma(2)}, \dots, i_{\sigma(n)})$$

The current implementation does **not** normalize by the product of hook lengths — each term has coefficient $\pm 1$ only. For a fully antisymmetric tableau `[1, 1, …, 1]`, the standard alternating sum is obtained.

**Parameters**:

| Parameter | Type | Description |
|---|---|---|
| `ctx` | `&'a AtomArena<'a>` | The expression arena |
| `tensor_expr` | `Atom<'a>` | The tensor expression (a `Fun` node) |
| `tableau` | `&YoungTableau` | The Young tableau shape |

**Returns**: `Atom<'a>` — the projected expression (an `Add` node whose terms carry signs).

**Example**:

```rust
use ocas_atom::tensor::young::{YoungTableau, young_project};
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// fully antisymmetric [1, 1]: f(a, b) → f(a, b) - f(b, a)
let f_ab = ctx.fun("f", &[ctx.var("a"), ctx.var("b")]);
let result = young_project(&ctx, f_ab, &YoungTableau::new(vec![1, 1]));
println!("{}", result);
// Output: f(a,b) - f(b,a)

// fully symmetric [2]: g(a, b) → g(a, b) + g(b, a)
let g_ab = ctx.fun("g", &[ctx.var("a"), ctx.var("b")]);
let result = young_project(&ctx, g_ab, &YoungTableau::new(vec![2]));
println!("{}", result);
// Output: g(a,b) + g(b,a)

// fully antisymmetric [1, 1, 1]: h(a, b, c) → alternating sum of 6 terms (coefficients ±1, not normalized)
let h = ctx.fun("h", &[ctx.var("a"), ctx.var("b"), ctx.var("c")]);
let result = young_project(&ctx, h, &YoungTableau::new(vec![1, 1, 1]));
println!("{}", result);
// Output: h(a,b,c) - h(b,a,c) - ... (6 terms, each with coefficient ±1)
```

**See also**: [YoungTableau](#youngtableau), [SymmetrySpec](#symmetryspec)

---

## Design Invariants

### `Atom` is a Copy arena handle

`Atom<'a>` is a Copy handle pointing to a node in the arena. `IndexSlot<'a>` is also Copy. This means tensors can be copied and compared cheaply — structural equality is equivalent to pointer equality (hash-consing).

### Explicit index matching

The tensor algebra of oCAS uses **explicit index matching** rather than the Einstein summation convention. Contractions must be performed manually by calling `contract`; indices are paired by label and position (upper/lower).

### Symmetry is advisory

The `Symmetry` enum and `SymmetrySpec` are metadata and are not used automatically by basic operations such as `contract`. Full symmetry handling requires:
1. Canonicalization (`canonicalize_tensors`) — symmetry enforcement at the level of graph isomorphism
2. Young projection (`young_project`) — explicit permutation expansion

### Graph isomorphism engine

Canonicalization in 0.22.0 is based on the McKay refinement-individualization graph isomorphism algorithm (implemented in the `tensor::graph` module). This engine is a standalone nauty implementation that provides:
- 1-WL color refinement (neighbor signatures + edge data + direction)
- equitable partitions + individualization-refinement DFS
- path-invariant pruning + automorphism orbits
- canonical form = the lexicographically maximal certificate

---

## Complete Example

```rust
use ocas_atom::tensor::{
    IndexPosition, IndexSlot, Symmetry, Tensor,
    contract, Contracted, symmetrise_sign,
};
use ocas_atom::tensor::spec::{SymmetrySpec, TensorRegistry};
use ocas_atom::tensor::canon::canonicalize_tensors;
use ocas_atom::tensor::dummy::refresh_dummies;
use ocas_atom::tensor::young::{YoungTableau, young_project};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);

    // 1. Create the metric tensor g_μν with symmetry
    let g = Tensor::new(Symbol::new("g"), vec![
        IndexSlot::new(ctx.var("mu"), IndexPosition::Lower),
        IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
    ]).with_symmetry(Symmetry::Symmetric);
    println!("g rank = {}, symmetry = {:?}", g.rank(), g.symmetry());
    // Output: g rank = 2, symmetry = Symmetric

    // 2. Contract g^μν A_ν → partial contraction
    let g_inv = Tensor::new(Symbol::new("g"), vec![
        IndexSlot::new(ctx.var("mu"), IndexPosition::Upper),
        IndexSlot::new(ctx.var("nu"), IndexPosition::Upper),
    ]);
    let a = Tensor::new(Symbol::new("A"), vec![
        IndexSlot::new(ctx.var("nu"), IndexPosition::Lower),
    ]);
    match contract(&ctx, &g_inv, &a) {
        Contracted::Product(tp) => {
            println!("Free indices after contraction: {}", tp.factors[0].rank());
            // Output: number of free indices after contraction: 1 (mu)
        }
        _ => {}
    }

    // 3. Canonicalize: j in T(i,j)*U(j,k) is renamed
    let mut reg = TensorRegistry::new();
    reg.register(Symbol::new("T"), SymmetrySpec::none());
    reg.register(Symbol::new("U"), SymmetrySpec::none());

    let prod = ctx.mul(&[
        ctx.fun("T", &[ctx.var("i"), ctx.var("j")]),
        ctx.fun("U", &[ctx.var("j"), ctx.var("k")]),
    ]);
    let ct = canonicalize_tensors(&ctx, prod, &reg).unwrap();
    println!("Canonical form: {}", ct.canonical_form);
    // Output: canonical form contains dummy indices d0, d1, etc.

    // 4. Young projection: fully antisymmetrize f(a, b)
    let f = ctx.fun("f", &[ctx.var("a"), ctx.var("b")]);
    let anti = young_project(&ctx, f, &YoungTableau::new(vec![1, 1]));
    println!("Antisymmetrized: {}", anti);
    // Output: f(a,b) - f(b,a)
}
```

---

## See also

- [Expression System](./rust-expressions.md) — the `Atom`, `AtomArena`, `Symbol` basic types
- [Rewriting and Simplification](./rust-rewrite.md) — pattern-matching-based expression transformation
- [Linear Algebra](../math/linear-algebra.md) — matrix operations (a low-rank special case of tensors)
- [Tensor Algebra and Canonicalization](../math/tensor-canonicalization.md) — the mathematical foundations
