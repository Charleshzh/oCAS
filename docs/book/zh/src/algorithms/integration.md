# 符号积分

oCAS 通过分层管线进行积分：快速启发式表、有理函数积分器、初等塔上
的 Risch 算法、三角到指数的重写，最后是特殊函数表。第一个产生答案
的层获胜。本章解释每一层以及何时返回未求值形式 `Integral(expr, var)`。

---

## 管线

`integrate(expr, var)` 按顺序尝试：

1. **启发式表** —— 幂法则、线性参数的 `sin`/`cos`/`exp`/`log`、线性
   替换。快速且总是最先尝试。
2. **有理函数积分器** —— Hermite 约化加对数部分（对数导数恒等式、
   配方、Rothstein–Trager）。处理 `x` 的任意有理函数。
3. **Risch 算法** —— 由 `log` 和 `exp` 构建的初等超越塔。
4. **三角重写** —— `sin`/`cos`/`tan`/… 重写为 `exp(I·x)` 后由 Risch
   积分，再尽力转换回实数形式。
5. **特殊函数表** —— 具有 `erf`、`Ei`、`Si`、Ci、Fresnel `S`/`C` 等
   闭式的非初等积分。
6. **未求值形式** —— `Integral(expr, var)`。

---

## 有理函数

积分变量的任意有理函数都被精确积分。多项式部分逐项积分；真分式由
Hermite 约化分解为有理部分加分母无平方的剩余；剩余部分产生对数
（经恒等式 `c·f'/f → c·log(f)`）、反正切（二次分母配方）或
Rothstein–Trager 对数。

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = parse(&ctx, "(2*x + 3)/(x^2 + 3*x + 5)").unwrap();
let result = integrate(&ctx, expr, Symbol::new("x"));
// log(x^2 + 3*x + 5)
```

---

## Risch 算法

初等超越被积函数通过构建*微分域塔* `ℚ(x, t₁, …, tₙ)` 处理，其中
每个 `tᵢ` 是下层域上的 `log` 或 `exp`，然后递归积分（Bronstein
《Symbolic Integration I》第 5 章）：

- 每层由 Hermite 约化分出有理部分；
- 对数部分使用对数导数恒等式；
- 多项式部分在 `log` 层用待定系数、在 `exp` 层用 Risch 微分方程
  `Dq + f·q = g` 积分；
- 基域 `ℚ(x)` 委托给有理函数积分器。

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
// ∫ x·exp(x) dx = (x - 1)·exp(x)
let result = integrate(&ctx, parse(&ctx, "x*exp(x)").unwrap(), Symbol::new("x"));
```

### 范围限制

当前片段只求 Risch 微分方程的**多项式**解，对数部分只使用对数导数
恒等式。因此：

- `∫ exp(x)/x dx` 没有初等原函数 —— 由特殊函数表回答为 `Ei(x)`。
- 某些需要自由选择常数使下层可积的 `log` 塔情形（如 `log(x+1)`）
  尚未判定，走回退。

当所有层都失败时，结果是未求值形式 `Integral(expr, var)`——这是
有意的答案，而非错误。

---

## 三角被积函数

`sin`、`cos`、`tan`、`cot`、`sec`、`csc` 经 `t = exp(I·x)` 重写为复
指数后由 Risch 积分。虚数单位作为常数塔生成元携带（`D I = 0`）。结果
在可能时转换回实数形式：共轭对数对合并为实 `log`/`atan` 项。

Risch 微分方程求解器目前在 `ℚ[x]` 上工作，因此系数含 `I` 的超指数
方程（如 `sin(x)·cos(x)` 或 `cos(x)²` 产生的方程）尚不能求解；这些被积
函数返回未求值形式。线性参数的简单 `sin`/`cos` 由启发式表覆盖。

---

## 特殊函数

没有初等原函数但有标准闭式的积分直接回答（定义与 SymPy 一致）：

| 被积函数 | 结果 |
|---|---|
| `exp(-x²)` | `(√π/2)·erf(x)` |
| `exp(x²)` | `(√π/2)·erfi(x)` |
| `exp(c·x²)`，`c < 0` | `√π/(2√(-c))·erf(√(-c)·x)` |
| `exp(x)/x` | `Ei(x)` |
| `sin(x)/x` | `Si(x)` |
| `cos(x)/x` | `Ci(x)` |
| `sinh(x)/x` | `Shi(x)` |
| `cosh(x)/x` | `Chi(x)` |
| `sin(x²)` | `√(π/2)·fresnels(√(2/π)·x)` |
| `cos(x²)` | `√(π/2)·fresnelc(√(2/π)·x)` |

```rust
use ocas::prelude::*;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
// ∫ exp(-x^2) dx = (√π/2)·erf(x)
let result = integrate(&ctx, parse(&ctx, "exp(-x^2)").unwrap(), Symbol::new("x"));
```

---

## 绑定

同一管线支撑 Python 与 C API：

- Python：`Expression.integrate(var)`
- C：`ocas_expr_integrate(...)`

找不到闭式时两者都返回未求值形式 `Integral(...)`，与 Rust API 一致。
