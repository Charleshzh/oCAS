# 求解器

oCAS 提供线性方程组、丢番图方程与多项式系统的求解器。本章介绍所有可用求解器及其用法。

---

## ℚ 上的线性方程组

`solve_linear_rational` 在有理数域上求解 $n \times n$ 系统 $Ax = b$。
输入系数为 `i64` 值；解以 `(分子, 分母)` 对返回。

```rust
let a = vec![vec![2, 1], vec![1, -1]];
let b = vec![5, 1];
let x = solve_linear_rational(&a, &b).unwrap();
// x = [(2, 1), (1, 1)]  → 2, 1
```

错误：`EmptySystem`、`NonSquare`、`Inconsistent`、`Underdetermined { rank }`。

Python：

```python
print(ocas.solve_linear_rational([[2, 1], [1, -1]], [5, 1]))
# [(2, 1), (1, 1)]
```

---

## ℤ 上的线性方程组

`solve_linear_integer` 求 $Ax = b$ 的整数解。若无整数解则返回错误。

```rust
// 2x + y = 3
let a = vec![vec![2, 1]];
let b = vec![3];
let x = solve_linear_integer(&a, &b).unwrap();
// x = [1, 1]  (2·1 + 1·1 = 3)
```

当解含分数时返回 `ResultNotInDomain` 错误。

---

## 丢番图方程

`solve_diophantine` 求解线性丢番图方程 $a \cdot x + b \cdot y = c$ 的整数解 $x, y$。

```rust
let sol = solve_diophantine(3, 5, 1).unwrap();
// sol = DiophantineSolution { x0: 2, y0: -1, x_step: 5, y_step: -3 }
```

结果给出特解 $(x_0, y_0)$ 和步长值。通解为：

$$
\begin{aligned}
x &= x_0 + x_{step} \cdot t \\
y &= y_0 + y_{step} \cdot t
\end{aligned}
$$

其中 $t$ 为任意整数。

---

## 多项式系统（基于 Gröbner 基）

`solve_polynomial_system` 首先计算 Gröbner 基，然后进行回代，求解多项式方程组。
它使用 Buchberger 算法，支持可配置的单项式序。

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// x + y = 0, x*y - 1 = 0  →  x + y = 0, y^2 + 1 = 0
let eq1 = parse(&ctx, "x + y").unwrap();
let eq2 = parse(&ctx, "x*y - 1").unwrap();
let sol = solve_polynomial_system(&ctx, &[eq1, eq2], &[Symbol::new("x"), Symbol::new("y")]);
```

结果为三角形多项式系统，可通过回代求解。

---

## 错误

所有求解器返回 `Result<T, SolveError>`。常见错误变体：

| 错误 | 含义 |
|---|---|
| `EmptySystem` | 未提供方程 |
| `NonLinear` | 系统对请求变量非线程 |
| `NonSquare` | 方程数与未知数个数不匹配 |
| `Inconsistent` | 无解 |
| `Underdetermined { rank }` | 无穷多解 |
| `ResultNotInDomain` | 解含分数但要求整数 |

---

## 参见

- [Rust API](./rust-api.md) — 系数域类型与多项式操作
- [重写与化简](./rewrite.md) — 化简求解结果
- [基准与性能对比](./performance.md) — Gröbner 基基准结果
