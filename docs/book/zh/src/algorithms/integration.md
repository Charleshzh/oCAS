# 符号积分

oCAS 通过分层管线进行积分：快速启发式表、有理函数积分器、初等塔上
的 Risch 算法、三角到指数的重写，最后是特殊函数表。第一个产生答案
的层获胜。本章解释每一层以及何时返回未求值形式 `Integral(expr, var)`。

---

## 管线

`integrate(expr, var)` 按顺序尝试：

1. **快速表** —— 结构内联的幂法则（$x^n$、$(ax+b)^n$，含分数指数）、
   线性参数 `sin`/`cos`/`exp` 的替换、以及 `log(x)` 等直接表项。快速且
   总是最先尝试。
2. **有理函数积分器** —— Hermite 约化加对数部分（对数导数恒等式、
   配方、Rothstein–Trager）。处理 `x` 的任意有理函数。
3. **符号常数有理积分器**（0.27）—— 系数在 ℚ(symbols) 中的有理函数：
   Yun 无平方分解 + Hermite + 部分分式（线性留数、二次因子的
   log/atan/atanh）。
4. **Risch 算法** —— 由 `log` 和 `exp` 构建的初等超越塔（塔深受
   `MAX_RISCH_DEPTH = 16` 上限约束）。
5. **三角重写** —— `sin`/`cos`/`tan`/… 重写为 `exp(I·x)` 后由 Risch
   积分，再尽力转换回实数形式。
6. **特殊函数表** —— 具有 `erf`、`Ei`、`Si`、Ci、Fresnel `S`/`C` 等
   闭式的非初等积分。
7. **规则表引擎**（0.27）—— 自研标准微积分规则（A–H 族：幂/二项式、
   指数/对数、三角、双曲、反三角/反双曲、有理拦截、根式、特殊形态），
   模板可含残项 `Integral(g, x)` 归约。可用
   `IntegrateOptions { rules: false }` 关闭。
8. **三角积化和差/降幂**（0.27.x）—— 线性变元 `sin`/`cos` 乘积与幂化为
   多倍角和式后逐项积分。
9. **启发式技巧** —— 分部积分、三角替换、Weierstrass $t = \tan(x/2)$、
   Euler 替换。
10. **有界展开重试**（0.27.x）—— 分配展开乘积为和式（项数预算 64）
    后逐项积分。
11. **未求值形式** —— `Integral(expr, var)`。

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

塔递归深度受 `MAX_RISCH_DEPTH = 16` 上限约束：超过该深度时 Risch 层
直接放弃，把被积函数交给下一层（防止病态被积函数造成无限递归）。

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

## Fuel 受限积分

`integrate_with_fuel` 包装同一管线，但将 [`Fuel`] 预算贯穿到两个积分后
化简阶段。病态结果（否则会导致重写器无限循环）可被确定性地截断。

```rust
use ocas_core::arena::Arena;
use ocas_core::fuel::Fuel;
use ocas_atom::{AtomArena, Symbol};
use ocas_calc::integral::integrate_with_fuel;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.var("x");
let fuel = Fuel::new(500);
let result = integrate_with_fuel(&ctx, expr, Symbol::new("x"), &fuel);
match result {
    Ok(r) => println!("{}", r),
    Err(_) => println!("fuel exhausted during simplification"),
}
```

仅在化简中途 fuel 耗尽时返回 `Err`。积分遍历本身使用内部深度限制；
fuel 仅约束化简后的传递。

---

## 绑定

同一管线支撑 Python 与 C API：

- Python：`Expression.integrate(var)`
- C：`ocas_expr_integrate(...)`

找不到闭式时两者都返回未求值形式 `Integral(...)`，与 Rust API 一致。
两者都**不认证**结果；需要可机检证书时用 `Expression.integrate_outcome`
（Python）或 `ocas_expr_integrate_outcome`（C）。

---

## 0.28.0：证书、循环检测与残项解析

### 符号证书与三值输出

`integrate_outcome`（Rust/Python/C 均已暴露）返回三值之一：

| 取值 | 含义 |
|---|---|
| `Found { value, certificate }` | 原函数 + **可机检的符号证书**（精确检查器中 `D(F) − f ≡ 0`） |
| `ProvedNonElementary { witness }` | 已证明非初等（0.28.0 无生产者，为后续非初等层预留） |
| `Unknown { residue, uncertified }` | 诚实未知：`residue` 是未求值 `Integral(f, x)`；`uncertified` 是管线候选（**不是答案**，仅供诊断） |

证书引擎（`ocas-calc/src/integral/certify.rs`）分层判定：结构零（`normalize` + 同类项收集）、
初等域零（三角→指数、`I² = −1` 约化、依赖 `exp`/`log` 原子合并后嵌入有理函数域）、根式域
（接口预留）。各层都是**可靠的**：把原子当作独立生成元只会漏判（假阴性），不会把错误答案判成
正确。域运算有确定性预算，超预算即诚实拒绝。

`ocas-calc/src/integral/risch.rs` 还会在**塔自身的域**里校验每个 Risch 结果，因此引擎不会发射
被精确检查器证伪的答案——0.28.0 正是靠这条护栏抓到并兜住一个潜伏错案（指数层有理部分对负幂次
系数的错误缩放）。

### 表达式级循环检测与残项解析

链入口（`integral/chain.rs`）按「表达式地址 + `rule_depth` + `parts_depth`」判环：**预算递减的
同形状重入**是流水线合法的递归下降，必须放行；只有预算相同的重复才是真循环。绝对条目兜底保持
0.27.3 的 256，残项解析的嵌套重入走独立预算，不与主链争用。

若干阶段会返回**部分结果**（分式分解已完成、剩下的 `Integral(...)` 本身可解）。顶层入口在链条
结束后会做一次预算化的残项解析：只处理「有理形状、≤2 个符号、≤64 节点」的残项，并且只在残项数
**严格减少**时接受解析结果。纯回退 `Integral(f, x)` 永远不解析（否则会重跑整条链，并在
`rules = false` 时违反调用方意图）。

### 依赖生成元合并（正则塔）

`ocas-calc/src/tower/merge.rs` 用精确域恒等式合并代数相关生成元：`exp(u)`/`exp(u+c)`、
`exp(u)`/`exp(−u)`、`log(u)`/`log(cu)`、`log(u)`/`log(u^k)`、`log(exp(u))`、`exp(log(u))`，
必要时登记常量生成元（`log(2)`、`exp(1)`，`D t = 0`）。塔现在返回**重写后的**被积表达式，
`risch` 在其中工作；判定不了的关系仍然诚实拒绝。
