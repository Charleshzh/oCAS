# Rust API

本章介绍 oCAS 的 Rust API。所有示例均假设：

```rust
use ocas::prelude::*;
```

在 `Cargo.toml` 中添加 `ocas = "0.13"`。

---

## 表达式

oCAS 的核心是由 arena 分配器管理的 `Atom` 表达式树。

```rust
let arena = Arena::new();
let ctx = AtomArena::new(&arena);

// 手动构建表达式
let x = ctx.var("x");
let y = ctx.var("y");
let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.mul(&[ctx.num(3), x]), ctx.num(1)]);

// 从字符串解析
let e = parse(&ctx, "sin(x)^2 + cos(x)^2").unwrap();
println!("{}", e);                       // sin(x)^2 + cos(x)^2
```

`Atom` 是可廉价复制的句柄（`Copy + Clone`）。arena 持有所有节点，在 drop 时一次性释放。

---

## 系数域

oCAS 通过 `Domain` 与 `EuclideanDomain` trait 支持多种系数域。

```rust
// 纯 Rust 大整数与有理数（默认构建）
let a = Integer::from(42);
let b = Integer::from(18);
let g = IntegerDomain.gcd(&a, &b);       // 6

// 有理数
let r = Rational::new(Integer::from(1), Integer::from(3));
println!("{}", r);                       // 1/3

// 有限域
let gf7 = FiniteField::new(7);
let fe = FiniteFieldElement::new(Integer::from(3), &gf7);
let inv = gf7.inv(&fe).unwrap();
println!("{}", inv);                     // 5  (3·5 ≡ 1 mod 7)

// 实数球算术（需 `mpfr` feature）
let ball = RealBallDomain.from_f64(1.0 / 3.0);
println!("{}", ball);                    // ~3.33333e-1 ± ε

// 复数
let z = Complex::new(RealBallDomain.from_f64(1.0), RealBallDomain.from_f64(2.0));
println!("{}", z);                       // (1.0 + 2.0i)
```

**`Assumptions`** 允许声明符号属性：

```rust
let mut assumptions = Assumptions::new();
assumptions.add("x", Assumption::Positive);
assumptions.add("n", Assumption::Integer);
assert!(assumptions.is_positive("x"));
```

---

## 多项式

### 稠密一元多项式

```rust
// x^2 + 3x + 2 在整数域
let p = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    &[Integer::from(2), Integer::from(3), Integer::from(1)],
);
println!("{}", p);                       // x^2 + 3*x + 2

// 求值
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

### 稀疏多元多项式

```rust
use std::collections::BTreeMap;

let mut terms = BTreeMap::new();
terms.insert(vec![1, 1], Integer::from(1));   // x*y
terms.insert(vec![2, 0], Integer::from(1));   // x^2
terms.insert(vec![0, 2], Integer::from(1));   // y^2
let sp = SparseMultivariatePolynomial::new(IntegerDomain, terms, Lex);
println!("{}", sp);                      // x^2 + x*y + y^2（lex 序）
```

可用单项式序：`Lex`（字典序）、`Grevlex`（分次反字典序）。

### 因式分解

```rust
// 无平方因子分解
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

对于整数域与有理数域上的完全因式分解，请使用 `factor` 模块（Hensel 提升 + 有限域方法）。

### Gröbner 基

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

### 根隔离

```rust
let p = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    &[Integer::from(-1), Integer::from(0), Integer::from(1)],  // x^2 - 1
);
for interval in p.real_root_intervals(&IntegerDomain) {
    println!("根在 [{}, {}]", interval.left(), interval.right());
}
```

---

## 矩阵

```rust
let m = Matrix::new(IntegerDomain, 2, 2, &[
    Integer::from(1), Integer::from(2),
    Integer::from(3), Integer::from(4),
]);

println!("{}", m.determinant());         // -2
println!("{}", m.rank());                // 2
println!("{}", m.trace());               // 5

// 转置
let mt = m.transpose();
assert_eq!(mt[(0, 1)], Integer::from(3));

// 矩阵乘法
let m2 = Matrix::new(IntegerDomain, 2, 2, &[
    Integer::from(2), Integer::from(0),
    Integer::from(0), Integer::from(2),
]);
let prod = m.matmul(&m2, &IntegerDomain);

// 在 ℚ 上解线性方程组
let a = Matrix::new(RationalDomain, 2, 2, &[
    Rational::from(2), Rational::from(1),
    Rational::from(1), Rational::from(1),
]);
let b = vec![Rational::from(4), Rational::from(3)];
let x = a.solve(&b, &RationalDomain).unwrap();
// x = [1, 2]
```

行列式使用 Bareiss 无分数算法。当行列式在域中可逆时，`inverse()` 返回精确逆。

---

## 微积分

```rust
let x = ctx.var("x");
let f = ctx.mul(&[ctx.num(2), ctx.pow(x, ctx.num(3))]);

// 微分
let df = diff(&ctx, f, Symbol::new("x"));
println!("{}", df);                      // 6*x^2

// Taylor 展开
let t = taylor(&ctx, f, Symbol::new("x"), ctx.num(0), 5);
println!("{}", t);                       // 2*x^3（对本多项式是精确的）

// 替换
let g = substitute(&ctx, f, x, ctx.add(&[x, ctx.num(1)]));
println!("{}", g);                       // 2*(x + 1)^3

// 积分（启发式）
let fi = integrate(&ctx, f, Symbol::new("x"));
println!("{}", fi);                      // 1/2*x^4
```

---

## 解析与输出

```rust
// 从字符串解析
let e = parse(&ctx, "x^2 + 2*x + 1").unwrap();

// Display 格式化输出中缀表达式
println!("{}", e);                       // x^2 + 2*x + 1

// 规范化（展平 Add/Mul、排序项、移除恒等元素）
let normalized = normalize(&ctx, e);
println!("{}", normalized);              // x^2 + 2*x + 1（已是规范形式）
```

解析器支持标准数学符号：`+`、`-`、`*`、`/`、`^`、括号、函数调用（`sin`、`cos`、`exp`、`log`、`sqrt`）及任意精度整数。

---

## 下一步

- [求解器](./solvers.md) — 线性方程组、丢番图方程、多项式系统
- [重写与化简](./rewrite.md) — 模式匹配、基于规则的化简
- [求值与 JIT](./evaluation.md) — 数值求值路径与性能
- [正确性](./correctness.md) — 交叉验证框架
