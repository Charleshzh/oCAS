# Automatic Differentiation

oCAS provides forward-mode automatic differentiation via *hyper-dual
numbers*. This is distinct from symbolic differentiation (`diff`): it
evaluates derivatives numerically alongside the function value, without
building or simplifying an expression tree.

---

## Hyper-Dual Numbers

A hyper-dual number extends the scalar field with infinitesimal
components ε₁, ε₂, … that satisfy εᵢ² = 0. Arithmetic on these
numbers propagates exact first (and optionally higher) partial
derivatives through every operation.

oCAS represents them as `HyperDual<T>` over any coefficient type `T`
that implements the `DualCoeff` trait. In practice `T` is typically
`Rational` (the only type with a complete `DualCoeff` implementation in
the standard build).

| Type | Role |
|---|---|
| `DualShape` | Layout descriptor: variable groupings, multiplication table |
| `HyperDual<T>` | Concrete dual number with value + derivative slots |
| `DualCoeff` | Trait: `zero`, `one`, arithmetic ops for `T` |
| `new_first_order(n)` | Convenience: shape with `n` first-order components |

---

## Quick Start

```rust
use std::sync::Arc;
use ocas_domain::Rational;
use ocas_domain::dual::{new_first_order, HyperDual};

// First-order dual for 2 variables
let shape = new_first_order(2);

// x = 1 + ε₁ (variable 0, unit coefficient)
let x = HyperDual::variable(&shape, 0, Rational::from(1));
// y = 2 + ε₂ (variable 1, unit coefficient)
let y = HyperDual::variable(&shape, 1, Rational::from(2));

// f = x * y  →  value = 2, ∂f/∂x = 2, ∂f/∂y = 1
let f = x * y;
println!("f    = {}", f.value());           // 2
println!("df/dx = {}", f.deriv(0).unwrap()); // 2
println!("df/dy = {}", f.deriv(1).unwrap()); // 1
```

---

## Available Operations

`HyperDual<T>` supports the standard arithmetic traits (`Add`, `Sub`,
`Mul`, `Div`, `Neg`). Each operation propagates derivatives via the
precomputed multiplication table in `DualShape`.

| Operation | Supported? |
|---|---|
| `+`, `-`, `*`, `/` | Yes |
| `inv()` (multiplicative inverse) | Yes |
| `pow` / `exp` / `log` / trig | **No** (transcendental functions require `DualCoeff` extensions not yet implemented) |

For transcendental functions, fall back to symbolic differentiation
(`diff`) or compile the expression to a numeric evaluator with
auto-diff support (planned post-1.0).

---

## Higher-Order Derivatives

`DualShape` supports multi-index groupings for higher-order derivatives
(e.g., second-order εᵢεⱼ components). Use `DualShape::new` with nested
index vectors for custom shapes, or the first-order convenience function
for the common case.

---

## Python & C Usage

### Python

```python
from ocas import DualShape, HyperDual, Rational

shape = DualShape.first_order(2)
x = HyperDual.variable(shape, 0, Rational(1))
y = HyperDual.variable(shape, 1, Rational(2))

f = x * y
print(f.value())    # 2
print(f.deriv(0))   # 2  (∂f/∂x)
print(f.deriv(1))   # 1  (∂f/∂y)

# Arithmetic dunder methods: __add__, __sub__, __mul__, __truediv__, __neg__
g = x + y * x
print(g.deriv(0))   # 3  (1 + y = 1 + 2)
```

### C

```c
#include <ocas.h>

ocas_OcasDualShape* shape = ocas_dual_shape_new(2, &err);
ocas_OcasHyperDual* x = ocas_dual_variable(shape, 0, "1", &err);
ocas_OcasHyperDual* y = ocas_dual_variable(shape, 1, "2", &err);

ocas_OcasHyperDual* f = ocas_dual_mul(x, y, &err);
char* val = ocas_dual_value(f, &err);    /* "2" */
char* dx  = ocas_dual_deriv(f, 0, &err); /* "2" */

ocas_string_free(val);
ocas_string_free(dx);
ocas_hyperdual_free(f);
ocas_hyperdual_free(y);
ocas_hyperdual_free(x);
ocas_dual_shape_free(shape);
```

See the [Python API](./bindings-python.md) and [C/C++ API](./bindings-c.md)
chapters for full documentation.

---

## Limitations

- Only `Rational` coefficients are supported in the standard build
  (no floating-point or finite-field duals yet).
- Transcendental functions (`sin`, `cos`, `exp`, `log`, `pow`) are not
  yet available on `HyperDual`.
- JIT compilation of dual expressions is not yet integrated.
