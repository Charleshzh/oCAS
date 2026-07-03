# Rust API

This chapter covers the Rust API of oCAS. All examples assume:

```rust
use ocas::prelude::*;
```

Import `ocas = "0.11"` in your `Cargo.toml`.

---

## Expressions

The core of oCAS is the `Atom` expression tree, managed by an arena allocator.

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// Building expressions manually
let x = ctx.var("x");
let y = ctx.var("y");
let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.mul(&[ctx.num(3), x]), ctx.num(1)]);

// Parsing from strings
let e = parse(&ctx, "sin(x)^2 + cos(x)^2").unwrap();
println!("{}", e);                       // sin(x)^2 + cos(x)^2
```

`Atom` is a cheaply-copyable handle (`Copy + Clone`). The arena owns all nodes
and deallocates them in one pass when dropped.

---

## Domains

oCAS supports multiple coefficient domains via the `Domain` and
`EuclideanDomain` traits.

```rust
// Pure-Rust big integers and rationals (default build)
let a = Integer::from(42);
let b = Integer::from(18);
let g = IntegerDomain.gcd(&a, &b);       // 6

// Rational numbers
let r = Rational::new(Integer::from(1), Integer::from(3));
println!("{}", r);                       // 1/3

// Finite fields
let gf7 = FiniteField::new(7);
let fe = FiniteFieldElement::new(Integer::from(3), &gf7);
let inv = gf7.inv(&fe).unwrap();
println!("{}", inv);                     // 5  (3·5 ≡ 1 mod 7)

// Real ball arithmetic (requires `mpfr` feature)
let ball = RealBallDomain.from_f64(1.0 / 3.0);
println!("{}", ball);                    // ~3.33333e-1 ± ε

// Complex numbers
let z = Complex::new(RealBallDomain.from_f64(1.0), RealBallDomain.from_f64(2.0));
println!("{}", z);                       // (1.0 + 2.0i)
```

**`Assumptions`** let you declare properties on symbols:

```rust
let mut assumptions = Assumptions::new();
assumptions.add("x", Assumption::Positive);
assumptions.add("n", Assumption::Integer);
assert!(assumptions.is_positive("x"));
```

---

## Polynomials

### Dense univariate

```rust
// x^2 + 3x + 2 over the integers
let p = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    &[Integer::from(2), Integer::from(3), Integer::from(1)],
);
println!("{}", p);                       // x^2 + 3*x + 2

// Evaluation
let val = p.evaluate(&Integer::from(5), &IntegerDomain);
println!("{}", val);                     // 42

// GCD
let q = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    &[Integer::from(1), Integer::from(1)],  // x + 1
);
let g = p.gcd(&q, &IntegerDomain);
println!("{}", g);                       // x + 1
```

### Sparse multivariate

```rust
use std::collections::BTreeMap;

let mut terms = BTreeMap::new();
terms.insert(vec![1, 1], Integer::from(1));   // x*y
terms.insert(vec![2, 0], Integer::from(1));   // x^2
terms.insert(vec![0, 2], Integer::from(1));   // y^2
let sp = SparseMultivariatePolynomial::new(IntegerDomain, terms, Lex);
println!("{}", sp);                      // x^2 + x*y + y^2 (lex order)
```

Available monomial orders: `Lex` (lexicographic), `Grevlex` (graded reverse lex).

### Factorization

```rust
// Square-free factorization
let p = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    &[Integer::from(-4), Integer::from(0), Integer::from(1)],  // x^2 - 4
);
for (factor, multiplicity) in p.square_free_factorization(&IntegerDomain) {
    println!("({})^ {}", factor, multiplicity);
}
// (x - 2)^1
// (x + 2)^1
```

For full factorization over integers and rationals, use the `factor` module
(Hensel lifting + finite-field methods).

### Gröbner bases

```rust
let x = ctx.var("x");
let y = ctx.var("y");
let polys = vec![
    SparseMultivariatePolynomial::from_coeffs(IntegerDomain, /* ... */),
    // ...
];
let basis: GroebnerBasis<Integer> = buchberger(&polys, &IntegerDomain, Grevlex);
for p in basis.polynomials() {
    println!("{}", p);
}
```

### Root isolation

```rust
let p = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    &[Integer::from(-1), Integer::from(0), Integer::from(1)],  // x^2 - 1
);
for interval in p.real_root_intervals(&IntegerDomain) {
    println!("root in [{}, {}]", interval.left(), interval.right());
}
```

---

## Matrices

```rust
let m = Matrix::new(IntegerDomain, 2, 2, &[
    Integer::from(1), Integer::from(2),
    Integer::from(3), Integer::from(4),
]);

println!("{}", m.determinant());         // -2
println!("{}", m.rank());                // 2
println!("{}", m.trace());               // 5

// Transpose
let mt = m.transpose();
assert_eq!(mt[(0, 1)], Integer::from(3));

// Matrix multiplication
let m2 = Matrix::new(IntegerDomain, 2, 2, &[
    Integer::from(2), Integer::from(0),
    Integer::from(0), Integer::from(2),
]);
let prod = m.matmul(&m2, &IntegerDomain);

// Solve linear system over ℚ
let a = Matrix::new(RationalDomain, 2, 2, &[
    Rational::from(2), Rational::from(1),
    Rational::from(1), Rational::from(1),
]);
let b = vec![Rational::from(4), Rational::from(3)];
let x = a.solve(&b, &RationalDomain).unwrap();
// x = [1, 2]
```

The determinant uses Bareiss's fraction-free algorithm. `inverse()` returns
the exact inverse when the determinant is invertible in the domain.

---

## Calculus

```rust
let x = ctx.var("x");
let f = ctx.mul(&[ctx.num(2), ctx.pow(x, ctx.num(3))]);

// Differentiation
let df = diff(&ctx, f, Symbol::new("x"));
println!("{}", df);                      // 6*x^2

// Taylor expansion
let t = taylor(&ctx, f, Symbol::new("x"), ctx.num(0), 5);
println!("{}", t);                       // 2*x^3 (exact for this polynomial)

// Substitution
let g = substitute(&ctx, f, x, ctx.add(&[x, ctx.num(1)]));
println!("{}", g);                       // 2*(x + 1)^3

// Integration (heuristic)
let fi = integrate(&ctx, f, Symbol::new("x"));
println!("{}", fi);                      // 1/2*x^4
```

---

## Parsing & Printing

```rust
// Parsing from string
let e = parse(&ctx, "x^2 + 2*x + 1").unwrap();

// Display formatting produces infix notation
println!("{}", e);                       // x^2 + 2*x + 1

// Normalization (flattens Add/Mul, sorts terms, removes identities)
let normalized = normalize(&ctx, e);
println!("{}", normalized);              // x^2 + 2*x + 1 (already canonical)
```

The parser supports standard mathematical notation: `+`, `-`, `*`, `/`, `^`,
parentheses, function calls (`sin`, `cos`, `exp`, `log`, `sqrt`), and
arbitrary-precision integers.

---

## Next Steps

- [Solvers](./solvers.md) — linear systems, Diophantine equations, polynomial systems
- [Rewrite & Simplification](./rewrite.md) — pattern matching, rule-based simplification
- [Evaluation & JIT](./evaluation.md) — numeric evaluation paths and performance
- [Correctness](./correctness.md) — cross-validation framework
