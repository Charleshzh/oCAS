# Getting Started / 快速上手

## Installation / 安装

**English**

oCAS is available as a Rust crate, a Python package, and a C/C++ library.

**中文**

oCAS 提供 Rust crate、Python 包和 C/C++ 库三种形式。

### Rust

Add to your `Cargo.toml`:

```toml
[dependencies]
ocas = "0.10"
```

```rust
use ocas::prelude::*;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let e = parse(&ctx, "x^2 + 2*x + 1").unwrap();
    let d = diff(&ctx, e, Symbol::new("x"));
    println!("{}", d); // 2*x + 2
}
```

### Python

```bash
pip install ocas
```

```python
import ocas

e = ocas.Expression("x^2 + 2*x + 1")
print(e.diff("x"))          # 2*x + 2
print(e.simplify())

# Polynomials over the integers
p = ocas.Polynomial([1, 2, 1])    # 1 + 2x + x^2
print(p.degree())                  # 2
print(p.eval(2))                   # "9"

# Matrices
m = ocas.Matrix([[1, 2], [3, 4]])
print(m.determinant())             # "-2"
print((m @ m).rows())

# Finite fields
gf5 = ocas.FiniteField(5)
q = ocas.Polynomial([1, 2, 1], domain=gf5)
print(q.eval(3))                   # "4"  (1 + 6 + 9 = 16 ≡ 4 mod 5)
```

### C/C++

Build the C library and link against `libocas_c`. See the
[C/C++ API](./bindings-c.md) chapter for details.

```c
#include <ocas.h>

ocas_expr* e = ocas_expr_parse("x^2 + 2*x + 1", NULL);
ocas_expr* d = ocas_expr_diff(e, "x", NULL);
char* s = ocas_expr_to_string(d, NULL);
printf("%s\n", s);   /* 2*x + 2 */
ocas_string_free(s);
ocas_expr_free(d);
ocas_expr_free(e);
```

---

## First Steps / 第一步

**English**

The most common entry points are:

| Task | Rust | Python |
|---|---|---|
| Parse an expression | `parse(&ctx, "x+1")` | `ocas.Expression("x+1")` |
| Differentiate | `diff(&ctx, e, x)` | `e.diff("x")` |
| Simplify | `simplify(&ctx, e, &rules, n)` | `e.simplify()` |
| Solve linear system | `solve_linear_rational(&a, &b)` | `ocas.solve_linear_rational(a, b)` |
| Numeric evaluation | `ExpressionEvaluator::<f64>` | `ocas.ExpressionEvaluator` |

**中文**

最常用的入口：

| 任务 | Rust | Python |
|---|---|---|
| 解析表达式 | `parse(&ctx, "x+1")` | `ocas.Expression("x+1")` |
| 微分 | `diff(&ctx, e, x)` | `e.diff("x")` |
| 化简 | `simplify(&ctx, e, &rules, n)` | `e.simplify()` |
| 解线性方程组 | `solve_linear_rational(&a, &b)` | `ocas.solve_linear_rational(a, b)` |
| 数值求值 | `ExpressionEvaluator::<f64>` | `ocas.ExpressionEvaluator` |
