# Tensors

oCAS provides basic tensor algebra with explicit index management. A
`Tensor` carries named index slots (upper/lower), optional symmetry
metadata, and can be converted to an `Atom` expression for symbolic
processing.

---

## Index Slots

Each index on a tensor has:

| Property | Type | Meaning |
|---|---|---|
| `label` | `Atom` | Index name (typically a Greek letter or `i`, `j`, …) |
| `position` | `IndexPosition` | `Upper` (contravariant) or `Lower` (covariant) |

```rust
use ocas_atom::tensor::{IndexSlot, IndexPosition};

let mu = IndexSlot::new("mu", IndexPosition::Upper);
let nu = IndexSlot::new("nu", IndexPosition::Lower);
```

---

## Creating Tensors

A `Tensor` is defined by a name (`Symbol`) and a vector of index slots:

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, Symmetry};

let slots = vec![
    IndexSlot::new("mu", IndexPosition::Upper),
    IndexSlot::new("nu", IndexPosition::Lower),
];
let t = Tensor::new(Symbol::new("g"), slots);

println!("rank = {}", t.rank());       // 2
println!("name = {}", t.name());       // g
```

Optionally attach symmetry metadata:

```rust
let symmetric = t.clone().with_symmetry(Symmetry::Symmetric);
let antisymmetric = t.clone().with_symmetry(Symmetry::Antisymmetric);
```

Symmetry metadata is advisory for downstream consumers; the tensor
algebra operations do not automatically symmetrize.

---

## Contraction

`contract` takes two tensors and returns either a scalar (if all indices
contract) or a `TensorProduct` with the remaining free indices:

```rust
use ocas_atom::tensor::{contract, Contracted};

// Contract two rank-1 tensors: A^μ B_μ → scalar
let a = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new("mu", IndexPosition::Upper),
]);
let b = Tensor::new(Symbol::new("B"), vec![
    IndexSlot::new("mu", IndexPosition::Lower),
]);

match contract(&ctx, &a, &b) {
    Contracted::Scalar(expr) => println!("scalar: {}", expr),
    Contracted::Product(tp) => {
        for f in &tp.factors {
            println!("free factor: {}", f.name());
        }
    }
}
```

Contraction matches indices by label and opposite position. Multiple
pairs of matching indices are contracted in a single call.

---

## Symmetrization

`symmetrise_sign` returns the sign of the antisymmetrization of a
tensor: `+1` for even permutations, `-1` for odd permutations of the
dummy (contracted) labels. This is useful for computing determinant-like
expressions.

```rust
use ocas_atom::tensor::symmetrise_sign;

let sign = symmetrise_sign(&tensor);
```

---

## Converting to Atoms

A tensor can be lowered to a standard `Atom` expression via `to_atom`:

```rust
let atom = tensor.to_atom(&ctx);
println!("{}", atom);  // e.g. g(mu, nu)
```

This allows tensors to participate in the general symbolic pipeline
(rewriting, simplification, differentiation).

---

## Python & C Usage

### Python

```python
from ocas import Tensor, contract_tensors, tensor_symmetrise_sign

# Create a rank-2 tensor
g = Tensor("g", [("mu", "upper"), ("nu", "lower")])
print(g.rank())    # 2
print(g.name())    # g

# Contract two tensors
A = Tensor("A", [("mu", "upper")])
B = Tensor("B", [("mu", "lower")])
result = contract_tensors(A, B)
print(result)      # scalar expression

# Antisymmetrization sign
sign = tensor_symmetrise_sign(g)
```

### C

```c
#include <ocas.h>

/* Create tensor A^μ */
const char* labels[] = {"mu"};
int positions[] = {1};  /* 1 = upper */
ocas_OcasTensor* A = ocas_tensor_create("A", labels, positions, 1, &err);

/* Query rank */
int rank = ocas_tensor_rank(A, &err);  /* 1 */

ocas_tensor_free(A);
```

See the [Python API](./bindings-python.md) and [C/C++ API](./bindings-c.md)
chapters for full documentation.

---

## Limitations

- Tensor algebra is explicit: indices are matched by label, not by
  Einstein summation convention. You must call `contract` manually.
- Graph-based canonicalization (automatic symmetry enforcement, Wick
  contractions) is deferred to post-1.0.
- Only `Symmetric` and `Antisymmetric` symmetry types are supported;
  no mixed or cyclic symmetries.
