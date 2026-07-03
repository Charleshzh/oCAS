# Solvers

oCAS provides solvers for linear systems, Diophantine equations, and
polynomial systems. This chapter covers the available solvers and their usage.

---

## Linear systems over ℚ

`solve_linear_rational` solves an $n \times n$ system $Ax = b$ over the
rational numbers. Input coefficients are `i64` values; the solution is
returned as `(numerator, denominator)` pairs.

```rust
let a = vec![vec![2, 1], vec![1, -1]];
let b = vec![5, 1];
let x = solve_linear_rational(&a, &b).unwrap();
// x = [(2, 1), (1, 1)]  → 2, 1
```

Errors: `EmptySystem`, `NonSquare`, `Inconsistent`, `Underdetermined { rank }`.

Python:

```python
print(ocas.solve_linear_rational([[2, 1], [1, -1]], [5, 1]))
# [(2, 1), (1, 1)]
```

---

## Linear systems over ℤ

`solve_linear_integer` finds integer solutions to $Ax = b$. It returns an
error if no integer solution exists.

```rust
// 2x + y = 3
let a = vec![vec![2, 1]];
let b = vec![3];
let x = solve_linear_integer(&a, &b).unwrap();
// x = [1, 1]  (2·1 + 1·1 = 3)
```

Errors include `ResultNotInDomain` when the solution involves fractions.

---

## Diophantine equations

`solve_diophantine` solves the linear Diophantine equation
$a \cdot x + b \cdot y = c$ for integer $x, y$.

```rust
let sol = solve_diophantine(3, 5, 1).unwrap();
// sol = DiophantineSolution { x0: 2, y0: -1, x_step: 5, y_step: -3 }
```

The result gives a particular solution $(x_0, y_0)$ and step values.
The general solution is:

$$
\begin{aligned}
x &= x_0 + x_{step} \cdot t \\
y &= y_0 + y_{step} \cdot t
\end{aligned}
$$

for any integer $t$.

---

## Polynomial systems (via Gröbner bases)

`solve_polynomial_system` solves systems of polynomial equations by computing
a Gröbner basis and then performing back-substitution. It uses the Buchberger
algorithm with configurable monomial order.

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// x + y = 0, x*y - 1 = 0  →  x + y = 0, y^2 + 1 = 0
let eq1 = parse(&ctx, "x + y").unwrap();
let eq2 = parse(&ctx, "x*y - 1").unwrap();
let sol = solve_polynomial_system(&ctx, &[eq1, eq2], &[Symbol::new("x"), Symbol::new("y")]);
```

The result is a simplified polynomial system in triangular form, which can
be solved by back-substitution.

---

## Errors

All solvers return `Result<T, SolveError>`. Common error variants:

| Error | Meaning |
|---|---|
| `EmptySystem` | No equations provided |
| `NonLinear` | System is not linear in the requested variables |
| `NonSquare` | Number of equations ≠ number of unknowns |
| `Inconsistent` | No solution exists |
| `Underdetermined { rank }` | Infinitely many solutions |
| `ResultNotInDomain` | Solution contains fractions when integers are required |

---

## See also

- [Rust API](./rust-api.md) — domain types and polynomial operations
- [Rewrite & Simplification](./rewrite.md) — simplifying solved expressions
- [Performance](./performance.md) — Gröbner basis benchmark results
