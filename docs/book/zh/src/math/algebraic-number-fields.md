# 高阶：代数数域

代数数域是计算机代数中连接多项式代数与数论的核心结构。在 oCAS 中，代数数域 $\mathbb{Q}(\alpha)$ 不仅是系数域——允许在扩域上做多项式运算和因式分解——还是 Galois 域 $\mathrm{GF}(p^d)$ 的统一抽象。本章从极小多项式和域扩张的基本概念出发，系统讲解代数数域上的算术、范数算法和模方法，并展示这些理论如何映射到 oCAS 的实现中。

## 前提知识

阅读本章前，建议先学习：

- [基础：多项式代数](./polynomial-algebra.md) — 多项式环、次数、带余除法、不可约性
- [基础：有限域与模算术](./finite-fields.md) — $\mathbb{F}_p$ 构造、模逆、CRT
- [基础：线性代数](./linear-algebra.md) — 矩阵运算、行列式、结式
- [进阶：多项式 GCD 与因式分解](./poly-gcd-factoring.md) — 无平方分解、Hensel 提升、模 GCD

### 极小多项式

设 $K$ 是 $\mathbb{Q}$ 的扩域，$\alpha \in K$。若存在非零多项式 $f \in \mathbb{Q}[x]$ 使得 $f(\alpha) = 0$，则称 $\alpha$ 为**代数元**（algebraic element）。满足此条件的**首一**最低次多项式 $m(x)$ 称为 $\alpha$ 的**极小多项式**（minimal polynomial）。

**性质**：

1. 极小多项式在 $\mathbb{Q}[x]$ 中**不可约**。
2. 任何满足 $f(\alpha) = 0$ 的多项式 $f$ 都有 $m \mid f$。
3. 极小多项式在同构意义下唯一。

**例**：$\sqrt{2}$ 的极小多项式是 $x^2 - 2$，因为 $(\sqrt{2})^2 - 2 = 0$ 且 $x^2 - 2$ 在 $\mathbb{Q}$ 上不可约（Eisenstein 判据，$p = 2$）。

### 域扩张与扩张次数

给定域 $F \subseteq K$，称 $K$ 是 $F$ 的一个**域扩张**（field extension），记为 $K/F$。$K$ 作为 $F$-向量空间的维数称为**扩张次数**：

$$[K : F] = \dim_F K$$

当 $K = F(\alpha)$（由单个代数元生成）时，称为**单扩张**（simple extension）。关键定理：

$$[F(\alpha) : F] = \deg(m_\alpha)$$

其中 $m_\alpha$ 是 $\alpha$ 在 $F$ 上的极小多项式。因此 $\{1, \alpha, \alpha^2, \ldots, \alpha^{d-1}\}$（$d = [K:F]$）构成 $K$ 作为 $F$-向量空间的一组基。

### 域塔与传递公式

若有域塔 $F \subseteq K \subseteq L$，则扩张次数满足**传递公式**：

$$[L : F] = [L : K] \cdot [K : F]$$

**例**：设 $\alpha = \sqrt[3]{2}$，$\omega = e^{2\pi i/3}$。则 $\mathbb{Q}(\alpha, \omega) / \mathbb{Q}(\alpha)$ 的次数为 2（$\omega$ 满足 $x^2 + x + 1$），$\mathbb{Q}(\alpha) / \mathbb{Q}$ 的次数为 3，因此 $[\mathbb{Q}(\alpha, \omega) : \mathbb{Q}] = 6$。

---

## 基础概念

### $\mathbb{Q}(\alpha)$ 的元素表示

设 $m(x) = x^d + a_{d-1}x^{d-1} + \cdots + a_0$ 为 $\alpha$ 的极小多项式，$d = [K:\mathbb{Q}]$。由 $m(\alpha) = 0$ 可得

$$\alpha^d = -(a_{d-1}\alpha^{d-1} + \cdots + a_0)$$

因此 $\mathbb{Q}(\alpha)$ 中的**每个元素**都可以唯一表示为

$$c_0 + c_1 \alpha + c_2 \alpha^2 + \cdots + c_{d-1} \alpha^{d-1}, \quad c_i \in \mathbb{Q}$$

即 $\alpha$ 的次数 $< d$ 的多项式。这个表示是**唯一的**，因为两个不同的次数 $< d$ 的多项式之差是 $m(x)$ 的一个非零倍数（次数 $< d$ 的多项式不可能是 $m$ 的倍数，除非为零）。

**加法**和**数乘**按系数逐项运算。**乘法**先做多项式乘法（得到次数 $\le 2(d-1)$ 的多项式），再对 $m(\alpha)$ 取模——即反复用 $\alpha^d = -\sum_{i=0}^{d-1} a_i \alpha^i$ 替换，直到次数 $< d$。

**例**：$\mathbb{Q}(\sqrt{2})$，$m(x) = x^2 - 2$。元素为 $a + b\sqrt{2}$（$a, b \in \mathbb{Q}$）。乘法：

$$(1 + \sqrt{2})(3 - \sqrt{2}) = 3 + 3\sqrt{2} - \sqrt{2} - 2 = 1 + 2\sqrt{2}$$

### 求逆

元素 $a(\alpha) \in \mathbb{Q}(\alpha)$ 可逆当且仅当 $a(\alpha) \ne 0$。求逆的关键工具是**扩展 Euclid 算法**：对 $a(x)$ 和 $m(x)$ 在 $\mathbb{Q}[x]$ 上做扩展 Euclid，得到

$$s(x) \cdot a(x) + t(x) \cdot m(x) = g(x)$$

其中 $g = \gcd(a, m)$。由于 $m$ 不可约且 $a \not\equiv 0 \pmod{m}$，必有 $\deg(g) = 0$，即 $g$ 是非零常数。归一化后 $g = 1$，故

$$a(\alpha)^{-1} = s(\alpha) \bmod m(\alpha)$$

**例**：在 $\mathbb{Q}(\sqrt{2})$ 中求 $(1 + \sqrt{2})^{-1}$。

用扩展 Euclid 求 $s(x)(1+x) + t(x)(x^2-2) = 1$：

- $x^2 - 2 = (x - 1)(1 + x) + (-1)$
- 因此 $(-1) = (x^2 - 2) - (x-1)(1+x)$
- 归一化：$1 = (x-1)(1+x) - (x^2-2)$

所以 $s(x) = x - 1$，验证：$(1+\sqrt{2})(\sqrt{2}-1) = \sqrt{2}^2 - 1 = 1$。✓

---

## 核心理论

### $\mathrm{GF}(p^d)$ 的构造

代数数域的模 $p$ "投影"需要有限域 $\mathrm{GF}(p^d)$。构造方法与 $\mathbb{Q}(\alpha)$ 完全平行：

1. 取 $\mathbb{F}_p[x]$ 中一个 $d$ 次不可约多项式 $\bar{m}(x)$（通常取极小多项式 $m(x)$ 模 $p$ 的像）。
2. $\mathrm{GF}(p^d) = \mathbb{F}_p[x]/(\bar{m}(x))$。
3. 元素为 $c_0 + c_1\alpha + \cdots + c_{d-1}\alpha^{d-1}$，$c_i \in \mathbb{F}_p$。
4. 运算与 $\mathbb{Q}(\alpha)$ 相同，只是系数在 $\mathbb{F}_p$ 中取模。

**关键条件**：$m(x) \bmod p$ 必须仍是 $d$ 次不可约多项式。若 $p$ 整除 $m$ 的某个分母或使 $m \bmod p$ 退化（降次或可约），则此素数不可用于模方法。

**例**：$m(x) = x^2 + 1$（$\mathbb{Q}(i)$ 的极小多项式）。在 $\mathbb{F}_3$ 上 $x^2 + 1$ 不可约（因为 $0^2+1=1$，$1^2+1=2$，$2^2+1=2$，均非零），所以 $\mathrm{GF}(3^2) = \mathbb{F}_3[x]/(x^2+1)$，乘法群阶为 8。

### Trager 范数算法

设 $K = \mathbb{Q}(\alpha)$，$f \in K[x]$ 是无平方多项式。**Trager 算法**通过将 $f$ 的因式分解问题从 $K[x]$ "降到" $\mathbb{Q}[x]$ 来求解。核心工具是**范数映射**。

#### 范数的定义

对 $f \in K[x]$，其**范数**（norm）$N(f) \in \mathbb{Q}[x]$ 定义为：

$$N(f)(x) = \operatorname{Res}_\alpha\bigl(m(\alpha),\, f(x, \alpha)\bigr)$$

即视 $f(x, \alpha)$ 为 $\alpha$ 的多项式（系数在 $\mathbb{Q}[x]$ 中），取它与 $m(\alpha)$ 的结式（消去 $\alpha$）。

**关键性质**：若 $f$ 在 $K[x]$ 中不可约，则 $N(f)$ 在 $\mathbb{Q}[x]$ 中是不可约多项式的幂。更精确地，若 $f = g_1 g_2 \cdots g_r$ 是 $K[x]$ 中的不可约分解，则

$$N(f) = N(g_1) \cdot N(g_2) \cdots N(g_r)$$

且每个 $N(g_i)$ 在 $\mathbb{Q}[x]$ 中不可约。

#### 求值—插值计算范数

直接计算结式 $\operatorname{Res}_\alpha(m, f(x,\alpha))$ 对高次多项式代价昂贵。oCAS 使用**求值—插值**法：

1. 范数 $N(f)$ 作为 $x$ 的多项式，次数为 $\deg_x(f) \cdot [K:\mathbb{Q}]$。
2. 取 $n = \deg_x(f) \cdot d + 1$ 个有理点 $x_0, x_1, \ldots, x_{n-1}$（使用对称序列 $0, -1, 1, -2, 2, \ldots$ 保持算术小）。
3. 在每个 $x_j$ 处，计算 $f(x_j, \alpha) \in K$（Horner 求值），再计算标量结式 $\operatorname{Res}_\alpha(m, f(x_j, \alpha))$。
4. 用 Newton 差商插值恢复 $N(f)$。

这是精确的（有理算术无舍入），且避免了高次结式的直接计算。

#### Trager 平移

即使 $f$ 无平方，其范数 $N(f)$ 也可能有平方因子（例如 $f = x - \alpha$ 在 $\mathbb{Q}(\sqrt{2})$ 上的范数是 $x^2 - 2$，已无平方；但 $f = x^2 - 2$ 的范数是 $(x^2-2)^2$，有平方）。

**Trager 平移**：对 $f(x)$ 做替换 $x \mapsto x - s\alpha$（$s \ge 0$），得到 $g(x) = f(x - s\alpha)$。由于坏的 $s$ 值（使 $N(g)$ 有平方的 $s$）只有有限多个，总存在 $s$ 使 $N(g)$ 无平方。

**算法**：

```
for s = 0, 1, 2, ...:
    g(x) ← f(x − s·α)          // Horner 复合
    R(x) ← N(g)(x)              // 求值—插值范数
    if R(x) 无平方（模素数验证）:
        return (s, g, R)
```

oCAS 尝试最多 `MAX_TRAGER_SHIFTS = 100` 个平移值。无平方性用模方法验证：取若干素数 $p$，若 $R \bmod p$ 在 $\mathbb{F}_p[x]$ 中无平方，则 $R$ 在 $\mathbb{Q}[x]$ 中无平方（充分条件，非必要——但这恰好是我们需要的"接受"判据）。

#### 因式分解流程

给定无平方 $f \in K[x]$：

1. **找平移**：$(s, g, R) \leftarrow \text{norm\_with\_shift}(f)$。
2. **分解范数**：$R = r_1 \cdot r_2 \cdots r_k$（在 $\mathbb{Q}[x]$ 中不可约分解）。
3. **回拉到 $K$**：对每个 $r_i$，计算 $\gcd_K(g, r_i)$（嵌入 $r_i$ 到 $K[x]$），得到 $f$ 在 $K[x]$ 中的因子。
4. **撤消平移**：$x \mapsto x + s\alpha$，恢复原变量。

**正确性保证**：范数的因子 $r_i$ 对应 $g$ 的因子，而 GCD 在 $K[x]$ 中精确提取这些因子。

### 数域上的模 GCD（Encarnación 方法）

第 3 步中的 $\gcd_K$（两个 $K[x]$ 多项式的 GCD）是计算瓶颈。直接在 $K[x]$ 上做 Euclidean GCD 会导致系数（$\mathbb{Q}(\alpha)$ 的元素，即有理数的 $\alpha$-多项式）剧烈膨胀。

**模方法**（Encarnación）将问题投影到 $\mathrm{GF}(p^d)$：

```
输入：a, b ∈ K[x]（K = ℚ(α)）
输出：monic gcd(a, b) ∈ K[x]

1. CRT 状态：residues = None, modulus = 1
2. for p in primes_from(2):
3.     if m(x) mod p 可约或降次: 跳过        // m 在 F_p 上退化，该素数不可用
4.     GF ← GF(p^d) = F_p[x]/(m mod p)
5.     if p | lc(a) 或 p | lc(b): 跳过
6.     g_p ← monic(gcd(a mod GF, b mod GF))   // GF(p^d)[x] 上的 GCD
7.     if deg(g_p) > best_deg: 跳过            // 不幸素数
8.     if deg(g_p) < best_deg:
9.         best_deg ← deg(g_p)
10.        residues ← g_p 的系数表  // [x-次数][α-次数] 的 F_p 值
11.        modulus ← p
12.    else:
13.        residues ← CRT(residues, g_p)  // 逐系数合并
14.        modulus ← modulus × p
15.    
16.    // 有理重构 + 试除验证
17.    candidate ← rational_reconstruct(residues, modulus)
18.    if candidate | a 且 candidate | b:
19.        return candidate
20.
21. return dense_euclidean_gcd(a, b)  // 安全回退
```

#### 素数筛选

不是所有素数都适用于模方法。需要满足：

1. **$m(x) \bmod p$ 保持不可约且次数不变**：否则 $\mathbb{F}_p[\alpha]/(m) \ne \mathrm{GF}(p^d)$，GCD 语义错误。
2. **$p$ 不整除首项系数**：否则映射后多项式降次。
3. **$p$ 不整除任何系数的分母**：有理数到 $\mathbb{F}_p$ 的映射失败。

#### 不幸素数检测

若素数 $p$ 上的模 GCD 次数**高于**当前最优值，该素数是"不幸的"——$p$ 整除了某个真因子的系数。直接丢弃此结果。

若模 GCD 次数**低于**当前最优值，说明之前的素数都是不幸的，丢弃之前的积累，以新结果重新开始。

#### 有理重构

CRT 合并后的残差是大整数（模大合数）。对每个系数（$\alpha$-多项式的每个分量），用**有理重构**从整数对中恢复有理数：

给定整数 $a$ 和模数 $m$，找有理数 $n/d$ 使得 $a \equiv n/d \pmod{m}$，且 $|n|, |d| \le \sqrt{m/2}$。

Wang/扩展 Euclid 方法追踪 $(r_i, t_i)$ 序列直到 $|r_1|, |t_1| \le \sqrt{m/2}$，验证 $a \cdot d \equiv n \pmod{m}$。

#### 试除验证

有理重构的候选 GCD 必须**试除验证**：检查它是否同时整除原始输入 $a$ 和 $b$。只有通过验证的结果才返回。

oCAS 最多尝试 `MAX_ANF_GCD_PRIMES = 1000` 个素数。极少需要超过数十个。极端情况下回退到 $K[x] = \mathbb{Q}(\alpha)[x]$ 上的稠密 Euclidean GCD（虽然慢但保证正确）。

### 分圆域 $\mathbb{Q}(\zeta_n)$

$n$ 次**本原单位根** $\zeta_n = e^{2\pi i/n}$ 生成的数域 $\mathbb{Q}(\zeta_n)$ 称为**分圆域**（cyclotomic field）。

#### 极小多项式

$\zeta_n$ 的极小多项式是**分圆多项式**（cyclotomic polynomial）：

$$\Phi_n(x) = \prod_{\substack{1 \le k \le n \\ \gcd(k,n)=1}} \left(x - e^{2\pi i k/n}\right)$$

其次数为 $\varphi(n)$（Euler 函数）。$\Phi_n$ 在 $\mathbb{Q}[x]$ 中不可约。

**例**：

- $\Phi_1(x) = x - 1$
- $\Phi_2(x) = x + 1$
- $\Phi_3(x) = x^2 + x + 1$
- $\Phi_4(x) = x^2 + 1$
- $\Phi_6(x) = x^2 - x + 1$
- $\Phi_p(x) = x^{p-1} + x^{p-2} + \cdots + 1$（$p$ 为素数）

分圆多项式满足 $x^n - 1 = \prod_{d \mid n} \Phi_d(x)$，可用 Möbius 反演逐因子计算。

#### 域的性质

- $[\mathbb{Q}(\zeta_n) : \mathbb{Q}] = \varphi(n)$。
- $\mathbb{Q}(\zeta_n)$ 是 $\mathbb{Q}$ 的**Abel 扩张**（Galois 群为 $(\mathbb{Z}/n\mathbb{Z})^*$）。
- 整数环为 $\mathbb{Z}[\zeta_n]$，判别式与 $n$ 的素因子有关。

**例**：$\mathbb{Q}(\zeta_5)$ 的扩张次数为 $\varphi(5) = 4$，极小多项式 $\Phi_5(x) = x^4 + x^3 + x^2 + x + 1$。元素表示为 $a_0 + a_1\zeta + a_2\zeta^2 + a_3\zeta^3$（$a_i \in \mathbb{Q}$）。

在 oCAS 中构造分圆域只需指定分圆多项式：

```rust
use ocas_domain::{AlgebraicNumberField, Domain, Rational, RationalDomain};

// Q(ζ_5): 极小多项式 x^4 + x^3 + x^2 + x + 1
// 升幂排列：[1, 1, 1, 1, 1]
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![
        Rational::new(1, 1),
        Rational::new(1, 1),
        Rational::new(1, 1),
        Rational::new(1, 1),
        Rational::new(1, 1),
    ],
);
let zeta = field.alpha(); // ζ_5
// ζ_5^5 = 1（由极小多项式保证）
let z5 = field.pow(&zeta, 5);
assert_eq!(z5, field.one());
```

---

## 在 oCAS 中的实现

### 核心结构体

代数数域的实现在 `ocas-domain/src/algebraic.rs` 中。

#### `AlgebraicElement<E>`

域元素——极小多项式的剩余类，以唯一的次数 $< d$ 的多项式代表存储：

```rust
pub struct AlgebraicElement<E> {
    coeffs: Vec<E>,  // 升幂排列，尾部零已裁剪
}
```

- 空向量表示零元素。
- `PartialEq` 直接比较系数向量（语义相等 = 结构相等，因为代表唯一）。
- `Display` 从高次到低次输出：`(c_d)·α^d + … + (c_1)·α + c_0`（常数项不带括号；例如 `1 + 2α` 显示为 `(2)·α + 1`）。

```rust
// 访问系数
let e = field.element(vec![r(1, 1), r(2, 1)]); // 1 + 2α
assert_eq!(e.coeffs(), &[r(1, 1), r(2, 1)]);
```

#### `AlgebraicExtension<D>`

域扩张 $D[\alpha]/(m(\alpha))$：

```rust
pub struct AlgebraicExtension<D: Domain> {
    base: D,              // 基域（RationalDomain 或 FiniteField）
    min_poly: Vec<D::Element>,  // 极小多项式（升幂，首一）
}
```

**构造**：

```rust
// Q(√2): m(x) = x² − 2，升幂排列 [−2, 0, 1]
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
```

**不变量**（`debug_assert` 检查）：
- 极小多项式次数 $\ge 1$（至少 2 个系数）。
- 极小多项式首一（首项系数 = 1）。
- **不验证不可约性**：可约模上的环有零因子，`inv` 对非单位返回 `None`。

**关键方法**：

| 方法 | 说明 |
|---|---|
| `extension_degree()` | $[K:\mathbb{Q}] = \deg(m)$ |
| `min_poly()` | 极小多项式系数（升幂） |
| `base_domain()` | 基域 $D$ |
| `alpha()` | 生成元 $\alpha$ |
| `from_base(c)` | 基域元素嵌入 |
| `element(coeffs)` | 从系数构造（自动对 $m$ 取模） |

#### `AlgebraicNumberField`

类型别名：

```rust
pub type AlgebraicNumberField = AlgebraicExtension<RationalDomain>;
```

### 算术运算

`AlgebraicExtension<D>` 实现了 `Domain` 和 `EuclideanDomain` trait。

**加法/减法**：多项式逐系数加减。

**乘法**：多项式乘法后对极小多项式取模（`reduce`）。取模过程：

```
while v.len() > d:
    c ← v.pop()                    // 最高次系数
    offset ← v.len() − d
    for i in 0..d:
        v[offset + i] −= c · m[i] // 减去 c·α^d 的展开
```

这利用了 $m$ 的首一性（首项系数 = 1），不需要除法。

**求逆**：扩展 Euclid 算法。计算 $s(x) \cdot a(x) + t(x) \cdot m(x) = 1$，返回 $s(\alpha) \bmod m(\alpha)$。若 $\gcd(a, m) \ne 1$（$a$ 是零因子），返回 `None`。

**除法**：`div(a, b) = a * inv(b)`。`inv(b)` 为 `None` 时除法失败。

**EuclideanDomain**：域上的带余除法退化为精确除法（余数恒为零）。GCD 退化为：两者均为零时返回零，否则返回一。

### 数域上的因式分解（Trager 算法）

因式分解的实现在 `ocas-poly/src/factor/algebraic.rs` 中。

#### 整体流程

```
输入：f ∈ K[x]（K = ℚ(α)）
输出：[(g₁, e₁), (g₂, e₂), …] 使得 f = lc(f) · ∏ gᵢ^eᵢ

1. if f 的系数全是 ℚ 常数:
       // 快速路径：先在 ℚ 上分解，再逐因子在 K 上分裂
       f_Q ← 去掉 α 分量（实际先去分母并取本原部分）
       for (h, e) in factor_primitive(f_Q):
           for g in factor_square_free_anf(K, embed(h)):
               output (g, e)
       return
2. for (g, e) in square_free_anf(K, f):      // Yun 无平方分解
3.     for h in factor_square_free_anf(K, g): // Trager 因式分解
4.         output (h, e)
```

#### 无平方分解（`square_free_anf`）

在 $K[x]$ 上实现 Yun 算法，但使用 `gcd_anf`（模 GCD）代替通用的伪余式 GCD。这是因为 $\mathbb{Q}(\alpha)$ 元素的系数是有理数，伪余式 GCD 的系数膨胀在任意精度有理数上更加严重。

#### Trager 因式分解（`factor_square_free_anf`）

```rust
fn factor_square_free_anf(field, f):
    if deg(f) <= 1:
        return [monic(f)]
    (s, g, R) ← norm_with_shift(field, f)  // 找平移使范数无平方
    rational_factors ← factor_square_free_rationals(R)  // ℚ 上分解范数
    if len(rational_factors) <= 1:
        return [monic(f)]  // 范数不可约 → f 不可约
    remaining ← g
    out ← []
    back_shift ← s·α  // 平移的逆
    for r_i in rational_factors:
        r_i_K ← embed(r_i)  // 嵌入 K[x]
        h ← gcd_anf(field, r_i_K, remaining)  // 模 GCD
        remaining ← remaining / h
        out.append(monic(compose_linear(h, back_shift)))  // 撤消平移
    return out
```

#### 模 GCD（`gcd_anf`）

Encarnación 方法的完整实现。CRT 状态跟踪每个系数的 $\alpha$-分量：

```rust
// CRT 状态：residues[x-次数][α-次数] 是整数残差
let mut residues: Option<Vec<Vec<Integer>>> = None;
let mut modulus = Integer::from(1);
```

每步迭代：

1. 取下一个素数 $p$，检查 $m \bmod p$ 是否不可约（`min_poly_irreducible_mod`：映射系数到 $\mathbb{F}_p$，验证次数和不可约性）。
2. 构造 $\mathrm{GF}(p^d) = \mathbb{F}_p[\alpha]/(m \bmod p)$。
3. 映射 $a, b$ 到 $\mathrm{GF}(p^d)[x]$，计算 monic GCD。
4. 合并到 CRT 状态。
5. 尝试有理重构 + 试除验证。

**不幸素数处理**：比较模 GCD 次数与历史最优值。高次 → 丢弃；低次 → 重置。

**回退**：超过 `MAX_ANF_GCD_PRIMES = 1000` 个素数后，使用 $K[x]$ 上的稠密 Euclidean GCD。

### 从 $\mathbb{Q}$ 到 $\mathrm{GF}(p^d)$ 的映射

将有理数映射到 $\mathbb{F}_p$：

```rust
fn map_rational_fp(c: &Rational, fp: &FiniteField) -> Option<FiniteFieldElement> {
    let num = fp.element(c.numer().to_bigint());  // 分子 mod p
    let den = fp.element(c.denom().to_bigint());  // 分母 mod p
    fp.div(&num, &den)  // num · den^{-1} mod p，分母为零时返回 None
}
```

将 $\mathbb{Q}(\alpha)$ 元素映射到 $\mathrm{GF}(p^d)$：

```rust
fn map_element_gf(e: &AlgebraicElement<Rational>, gf: &GaloisField) -> Option<...> {
    e.coeffs().iter().map(|c| map_rational_fp(c, gf.base_domain())).collect()
}
```

映射 $\alpha$-多项式到 $\mathrm{GF}(p^d)[x]$：逐系数映射。

### 平移与复合

`compose_linear(f, shift)` 用 Horner 方法计算 $f(x + \text{shift})$：

```rust
fn compose_linear<D: Domain>(f: &UP<D>, shift: &D::Element) -> UP<D> {
    let domain = f.domain().clone();
    // linear = x + shift
    let linear = UP::from_coeffs(domain.clone(), vec![shift.clone(), domain.one()]);
    let mut acc = UP::from_coeffs(domain.clone(), Vec::new()); // 零多项式
    for c in f.coeffs().iter().rev() {
        // acc ← acc · (x + shift) + c（多项式 Horner 复合）
        acc = acc
            .mul(&linear)
            .add(&UP::from_coeffs(domain.clone(), vec![c.clone()]));
    }
    acc
}
```

在 Trager 算法中，平移为 $x \mapsto x - s\alpha$，逆平移为 $x \mapsto x + s\alpha$。

### 范数计算

`norm_eval_interp` 实现求值—插值范数：

```rust
fn norm_eval_interp(field: &AlgebraicNumberField, g: &UP<AlgebraicNumberField>) -> UP<RationalDomain> {
    let deg_x = g.degree().unwrap_or(0);
    let d = field.extension_degree();
    let n_points = d * deg_x + 1;  // 插值需要的点数
    // 对称点 0, −1, 1, −2, 2, ...
    for j in 0..n_points {
        let x = Rational::new(((j+1)/2) * sign, 1);
        // Horner 求值 g(x) ∈ K
        let val = horner_eval(g, field.from_base(x));
        // 标量结式 Res_α(m, val(α))
        let norm_val = m_poly.resultant(&val_as_poly);
        xs.push(x); ys.push(norm_val);
    }
    interpolate_rational(&xs, &ys)  // Newton 差商插值
}
```

### Galois 域的统一抽象

`AlgebraicExtension` 的泛型设计使得 $\mathbb{Q}(\alpha)$ 和 $\mathrm{GF}(p^d)$ 共享同一份代码：

| 类型 | 基域 $D$ | 极小多项式 | 用途 |
|---|---|---|---|
| `AlgebraicExtension<RationalDomain>` | $\mathbb{Q}$ | $\mathbb{Q}[x]$ 中不可约 | 精确算术 |
| `AlgebraicExtension<FiniteField>` | $\mathbb{F}_p$ | $\mathbb{F}_p[x]$ 中不可约 | 模方法的"工作域" |

在 `gcd_anf` 中，同一素数 $p$ 上的 GCD 计算在 $\mathrm{GF}(p^d)[x]$ 中完成，然后 CRT 合并回 $\mathbb{Q}(\alpha)[x]$。

---

## 进阶话题

### 可约模与零因子

当 $m(x) \bmod p$ 在 $\mathbb{F}_p$ 上可约时，$\mathbb{F}_p[x]/(m)$ 不是域——它有零因子。oCAS 不验证极小多项式的不可约性，允许在可约模上构造环。此时：

- 非零零因子的 `inv` 返回 `None`（扩展 Euclid 的 $\gcd \ne 1$）。
- `gcd_anf` 跳过此类素数（`min_poly_irreducible_mod` 检查）。

**例**：$m(x) = x^2 - 1 = (x-1)(x+1)$ 在 $\mathbb{Q}$ 上可约。$\mathbb{Q}[x]/(x^2-1)$ 有零因子 $\alpha - 1$（因为 $(\alpha-1)(\alpha+1) = \alpha^2-1 = 0$）。

```rust
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-1, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
let alpha = field.alpha();
let a = field.sub(&alpha, &field.one());
assert!(field.inv(&a).is_none()); // α − 1 是零因子
```

### 模 GCD 的素数筛选效率

在实践中，大多数小素数都满足 $m \bmod p$ 不可约的条件。对于度 $d$ 的极小多项式，$\mathbb{F}_p$ 上随机 $d$ 次多项式不可约的概率约为 $1/d$（由不可约多项式计数公式）。因此平均只需尝试约 $d$ 个素数就能找到一个可用的。

oCAS 从 $p = 2$ 开始遍历素数（使用 `primes_from` 迭代器）。对于较大的 $d$（如 $d \ge 6$），可能需要跳过较多素数，但 1000 个素数的预算几乎总是足够的。

### 系数膨胀与有理重构的权衡

模方法的核心优势是避免了系数膨胀：$\mathrm{GF}(p^d)$ 上的运算是模运算，系数有界。但 CRT 合并后的残差增长为 $O(\prod p_i)$，有理重构需要 $\sqrt{m/2}$ 的界。

实践中，若真 GCD 的系数较小（典型情况），通常 2–5 个素数就够了。若系数很大（如高次扩张上的复杂 GCD），可能需要更多素数来积累足够的 CRT 模数。

### 与多变量因式分解的关系

Trager 范数算法也可以用于**多变量**多项式在代数数域上的因式分解（Wang EEZ + Zassenhaus 重组），但这在 oCAS 中尚未实现。当前的多元因式分解（`ocas-poly/src/factor/multivariate.rs`）仅支持 $\mathbb{Z}$ 和 $\mathbb{F}_p$ 系数。

### 计算复杂度

对 $f \in K[x]$（$K = \mathbb{Q}(\alpha)$，$[K:\mathbb{Q}] = d$，$\deg f = n$）：

| 步骤 | 复杂度 | 说明 |
|---|---|---|
| 范数计算 | $O(dn \cdot (n + d^3))$ 次系数运算 | 每个求值点：Horner 求值 $O(n)$ + 结式 $O(d^3)$；共 $O(dn)$ 个点 |
| 无平方分解 | $O(n^2)$ 次 GCD 调用 | Yun 算法，瓶颈在 GCD |
| 模 GCD | $O(d \cdot p \cdot \log p)$ 每个素数 | $\mathrm{GF}(p^d)$ 上的 Euclidean GCD |
| CRT 合并 | $O(d \cdot n \cdot \log^2 M)$ | $M$ 为 CRT 模数，$n$ 个系数 × $d$ 个 $\alpha$-分量 |

总体瓶颈通常是范数计算中的结式求值（$d$ 较大时）或模 GCD 的素数筛选（不幸素数较多时）。

---

## 参考文献

- **[Cohen]** H. Cohen, *A Course in Computational Algebraic Number Theory*, Springer GTM 138, 1993. — 第 2 章（数域的表示与算术）、第 3 章（基本数论函数）、第 4 章（理想与类群）。
- **[Trager]** B. M. Trager, "Algebraic factoring and rational function integration", *Proc. SYMSAC 1976*, pp. 196–208. — Trager 范数算法的原始论文。
- **[Encarnación]** M. J. Encarnación, "Computing GCDs of polynomials over algebraic number fields", *J. Symbolic Computation* 20:3, 1995, pp. 299–313. — 数域上模 GCD 方法。
- **[Gathen–Gerhard]** J. von zur Gathen & J. Gerhard, *Modern Computer Algebra*, 3rd ed., Cambridge, 2013. — 第 6 章（域扩张）、第 14–15 章（因式分解）。
- **[Cox–Little–O'Shea]** D. Cox, J. Little & D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed., Springer UTM, 2015. — 第 4 章（多项式环上的 Groebner 基与消元）。
- **[Lang]** S. Lang, *Algebra*, 3rd ed., Springer GTM 211, 2002. — 第 5 章（域扩张的基本理论）。

**参见**：[系数域 API](../api/rust-domains.md) — `AlgebraicExtension`、`AlgebraicNumberField`、`AlgebraicElement` 完整 API；[因式分解 API](../api/rust-factoring.md) — `DenseUnivariatePolynomial::factor()` 在代数数域上的使用；[多项式 GCD 与因式分解](./poly-gcd-factoring.md) — 模 GCD 和 Hensel 提升的通用理论；[有限域与模算术](./finite-fields.md) — $\mathbb{F}_p$ 和 $\mathrm{GF}(p^d)$ 的构造。
