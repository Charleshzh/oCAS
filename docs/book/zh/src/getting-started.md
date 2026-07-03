# 快速上手

## 安装

oCAS 提供 Rust crate、Python 包和 C/C++ 库三种形式。

### Rust

在你的 `Cargo.toml` 中添加：

```toml
[dependencies]
ocas = "0.11"
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

# 整数多项式
p = ocas.Polynomial([1, 2, 1])    # 1 + 2x + x^2
print(p.degree())                  # 2
print(p.eval(2))                   # "9"

# 矩阵
m = ocas.Matrix([[1, 2], [3, 4]])
print(m.determinant())             # "-2"
print((m @ m).rows())

# 有限域
gf5 = ocas.FiniteField(5)
q = ocas.Polynomial([1, 2, 1], domain=gf5)
print(q.eval(3))                   # "4"  (1 + 6 + 9 = 16 ≡ 4 mod 5)
```

### C/C++

构建 C 库并链接 `libocas_c`。详见 [C/C++ API](./bindings-c.md) 章节。

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

## 第一步

最常用的入口：

| 任务 | Rust | Python |
|---|---|---|
| 解析表达式 | `parse(&ctx, "x+1")` | `ocas.Expression("x+1")` |
| 微分 | `diff(&ctx, e, x)` | `e.diff("x")` |
| 化简 | `simplify(&ctx, e, &rules, n)` | `e.simplify()` |
| 解线性方程组 | `solve_linear_rational(&a, &b)` | `ocas.solve_linear_rational(a, b)` |
| 数值求值 | `ExpressionEvaluator::<f64>` | `ocas.ExpressionEvaluator` |
