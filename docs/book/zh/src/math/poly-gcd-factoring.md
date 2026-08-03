# 进阶：多项式 GCD 与因式分解

本章系统讲解多项式最大公因式（GCD）计算和因式分解的核心算法。这些算法是计算机代数系统的基础构件——Gröbner 基计算、有理函数化简、ODE 求解中的部分分式分解都依赖于它们。我们从经典的 Euclidean 算法出发，逐步引入模方法、Hensel 提升和多元扩展，最终覆盖有限域和代数数域上的因式分解。

## 前提知识

阅读本章前，建议先学习：

- [基础：多项式代数](./polynomial-algebra.md) — 多项式环、次数、单项式序
- [基础：有限域与模算术](./finite-fields.md) — $\mathbb{F}_p$ 构造、模逆、CRT
- [基础：线性代数](./linear-algebra.md) — 矩阵运算、行列式

### Euclidean 算法回顾

设 $\mathbb{F}$ 为域，$f, g \in \mathbb{F}[x]$，$g \ne 0$。Euclidean 算法通过反复带余除法计算 GCD：

$$f = q_1 g + r_1, \quad g = q_2 r_1 + r_2, \quad r_1 = q_3 r_2 + r_3, \quad \ldots, \quad r_{k-2} = q_k r_{k-1} + r_k$$

最后一个非零余式 $r_k$ 就是 $\gcd(f, g)$（精确到常数倍）。

在 $\mathbb{Z}[x]$ 中，系数不再是域元素，带余除法可能不可行。我们需要**伪除法**（pseudo-division）：

$$\text{lc}(g)^{d+1} \cdot f = q \cdot g + r, \quad d = \deg(f) - \deg(g)$$

其中 $\text{lc}(g)$ 是 $g$ 的首项系数。这保证了商和余式都在 $\mathbb{Z}[x]$ 中，但代价是余式的系数可能膨胀。

### 伪余式序列（PRS）

直接对 $\mathbb{Z}[x]$ 使用伪除法的 Euclidean 算法会导致**系数膨胀**：在度约 16 以上时，中间结果的系数可能达到数百位。为控制膨胀，文献中提出了多种**伪余式序列**（Pseudo-Remainder Sequences）：

- **原始 PRS**（primitive PRS）：每步取伪余式的**本原部分**（除以系数的内容），使中间余式互素、避免不必要的系数膨胀
- **子结式 PRS**（Subresultant PRS）：通过精心选择的缩放因子，使中间系数的大小恰好等于子结式行列式的理论界，避免了原始 PRS 的过度膨胀

在 oCAS 中，`DenseUnivariatePolynomial::gcd` 使用带伪余式的 Euclidean 算法（`ocas-poly/src/gcd.rs`），返回本原 GCD。伪余式在度数较高时（实践中约 $\deg > 16$）系数急剧膨胀，此时应改用 `gcd_modular_z`（Brown 模 GCD，`ocas-poly/src/gcd/modular.rs`）；目前两者由调用方显式选择，没有自动的度数调度。

### 结式（Resultant）

两个多项式 $f, g \in \mathbb{F}[x]$ 的**结式** $\text{Res}(f, g)$ 是一个标量，满足：

$$\text{Res}(f, g) = 0 \iff \gcd(f, g) \ne 1$$

结式可以通过 Sylvester 矩阵的行列式计算，也可以通过 PRS 高效计算（Brown 的 PRS 算法）。在 oCAS 中，`resultant()` 方法使用 Brown PRS 算法实现。

## 基础概念

### GCD 的定义

设 $D$ 为唯一分解整环（UFD），$f, g \in D[x]$。$f$ 和 $g$ 的**最大公因式** $\gcd(f, g)$ 是满足以下条件的首一多项式 $d$：

1. $d \mid f$ 且 $d \mid g$（公因式）
2. 对任意公因式 $c$，有 $c \mid d$（最大性）

在 $\mathbb{F}[x]$（$\mathbb{F}$ 为域）中，GCD 总是首一的（乘以首项系数的逆）。在 $\mathbb{Z}[x]$ 中，GCD 通常取为**本原的**（primitive），即系数的 GCD 为 1。

### 内容与本原化

多项式 $f = \sum a_i x^i \in \mathbb{Z}[x]$ 的**内容**（content）定义为系数的 GCD：

$$\text{cont}(f) = \gcd(a_0, a_1, \ldots, a_n)$$

**本原部分**（primitive part）定义为：

$$\text{pp}(f) = f / \text{cont}(f)$$

**Gauss 引理**：两个本原多项式的乘积仍是本原的。由此可得：

$$\gcd(f, g) = \gcd(\text{cont}(f),\, \text{cont}(g)) \cdot \gcd(\text{pp}(f),\, \text{pp}(g))$$

这允许我们将内容（整数 GCD）和本原部分（多项式 GCD）分开计算。

在 oCAS 中：

```rust
use ocas_poly::DenseUnivariatePolynomial;
use ocas_domain::IntegerDomain;

let f = DenseUnivariatePolynomial::from_coeffs(
    IntegerDomain,
    vec![Integer::from(6), Integer::from(4), Integer::from(2)],
);
assert_eq!(f.content(), Integer::from(2));       // cont(f) = gcd(6,4,2) = 2
assert_eq!(f.primitive_part().coeffs(),           // pp(f) = 3 + 2x + x²
    &[Integer::from(3), Integer::from(2), Integer::from(1)]);
```

## 核心理论

### 无平方分解（Square-free Factorization）

#### 动机

完全因式分解的第一步是将多项式分解为无平方因子的幂之积：

$$f = c \cdot g_1^1 \cdot g_2^2 \cdot g_3^3 \cdots$$

其中 $g_i$ 互素且无平方因子（$\gcd(g_i, g_i') = 1$）。这大幅简化了后续的不可约分解——每个 $g_i$ 只需分解一次。

#### Yun 算法

对于特征 0 的域（如 $\mathbb{Q}$ 或 $\mathbb{Z}$），Yun 算法基于以下观察：若 $f = g_1 g_2^2 g_3^3 \cdots$，则

$$\gcd(f, f') = g_2 g_3^2 g_4^3 \cdots$$

**算法步骤**：

```
输入：f ∈ ℤ[x]（已本原化）
输出：[(g₁, 1), (g₂, 2), …] 使得 f = ∏ gₖᵏ

1. f ← pp(f),  f' ← derivative(f)
2. g ← gcd(f, f')
3. w ← f / g
4. k ← 1
5. while w ≠ 1:
6.     h ← gcd(w, g)
7.     z ← w / h          // z = gₖ（重数为 k 的无平方因子）
8.     if z ≠ 1: 输出 (z, k)
9.     g ← g / h
10.    w ← h
11.    k ← k + 1
```

**关键性质**：每步迭代的 GCD 计算精确消去了更高重数的因子，使得 $z = w/h$ 恰好包含重数 $k$ 的因子。

**复杂度**：$O(n^2)$ 次系数运算（$n = \deg f$），瓶颈在于 GCD 计算。

#### 有限域上的特殊情况

在 $\mathbb{F}_p[x]$ 中，当 $p \mid \deg f$ 时，形式导数 $f'$ 可能为零（例如 $f(x) = x^p - 1$ 在 $\mathbb{F}_p$ 上有 $f' = 0$）。此时 $f$ 是某个多项式的 $p$ 次幂：$f(x) = g(x^p) = g(x)^p$。Musser/Bernardin 算法检测此情况，取 $p$ 次根后递归，再将重数乘以 $p$。

在 oCAS 中，`DenseUnivariatePolynomial::square_free_factorization()` 实现了 Yun 算法（`ocas-poly/src/factor/mod.rs`），`square_free_factorization_ff()` 实现了有限域版本（`ocas-poly/src/factor/finite_field.rs`）。

### 模 GCD（Brown 1971）

#### 问题

直接在 $\mathbb{Z}[x]$ 上用伪除法计算 GCD 会导致系数指数级膨胀。Brown 的模算法将问题"投影"到多个 $\mathbb{F}_p[x]$ 上，利用域上的高效 GCD，再用 CRT 重建。

#### 算法核心思想

设 $a, b \in \mathbb{Z}[x]$，令 $g = \gcd(a, b)$（本原）。关键观察：

$$\text{lc}(g) \mid \gcd(\text{lc}(a), \text{lc}(b)) \equiv \gamma$$

因此 $\gamma^{-1} \cdot \text{lc}(g) \mid 1$，即模 $p$ 的 monic GCD 乘以 $\gamma$ 后给出 $\mathbb{Z}[x]$ 上 GCD 的一个缩放版本。

#### 算法步骤

```
输入：a, b ∈ ℤ[x]（非零）
输出：gcd(a, b) ∈ ℤ[x]（本原）

1. a_p ← pp(a),  b_p ← pp(b)
2. γ ← gcd(lc(a_p), lc(b_p))    // GCD 首项系数的倍数
3. images ← [],  best_deg ← ∞
4. for p in primes_from(2³⁰):    // 从大于 2³⁰ 的素数开始
5.     if p | γ: continue
6.     field ← FiniteField(p)
7.     g_p ← monic(gcd(a_p mod p, b_p mod p))
8.     g_scaled ← g_p · (γ mod p) · (lc(g_p)⁻¹ mod p)
9.     
10.    if deg(g_scaled) > best_deg: continue    // 不幸素数，丢弃
11.    if deg(g_scaled) < best_deg:             // 发现更小的 GCD
12.        best_deg ← deg(g_scaled)
13.        images ← [(p, g_scaled)]
14.    else: images.append((p, g_scaled))
15.    
16.    // CRT 重建 + 试除验证
17.    candidate ← primitive_part(CRT_reconstruct(images))
18.    if candidate divides a_p AND candidate divides b_p:
19.        return candidate
20.    
21. if 已尝试的素数个数 > MAX_PRIMES: break
22.
23. return fallback_pseudo_remainder_gcd(a, b)  // 极端回退
```

#### 关键细节

**素数选择**：素数从 $> 2^{30}$ 开始。每个素数在 CRT 重建中贡献约 30 比特的精度，因此系数位长为 $L$ 的多项式通常只需约 $L/30$ 个素数即可重建（实践中一般不超过几十个）。

**Landau-$\gamma$ 缩放**：monic 模 GCD $g_p$ 需要乘以 $\gamma = \gcd(\text{lc}(a), \text{lc}(b))$ 才能与 $\mathbb{Z}[x]$ 上的 GCD 对齐。具体地，$\text{lc}(g_p)^{-1} \cdot \gamma \bmod p$ 是缩放因子。

**CRT 对称代表**：重建时使用对称剩余 $(-p/2, p/2]$ 而非 $[0, p)$，使系数绝对值最小化。

**不幸素数检测**：若某素数 $p$ 上的模 GCD 次数高于当前最优值，该素数是"不幸的"（$p$ 整除真 GCD 的某个系数的分母），直接丢弃。若次数更低，说明之前的素数都是不幸的，清空重来。

**试除验证**：CRT 重建的候选 $c$ 必须满足 $c \mid a_p$ 且 $c \mid b_p$（精确除法），否则继续收集更多素数。

**安全回退**：尝试最多 10000 个素数后（实践中极少超过数十个），回退到伪余式 GCD 保证终止。

在 oCAS 中，`gcd_modular_z()` 实现在 `ocas-poly/src/gcd/modular.rs`：

```rust
use ocas_poly::gcd::modular::gcd_modular_z;

let a = DenseUnivariatePolynomial::from_coeffs(d, vec![i(-1), i(0), i(1)]); // x² - 1
let b = DenseUnivariatePolynomial::from_coeffs(d, vec![i(1), i(2), i(1)]);  // x² + 2x + 1
let g = gcd_modular_z(&a, &b);
assert_eq!(g.coeffs(), &[i(1), i(1)]);  // x + 1
```

### 多元 GCD（Brown 求值—插值）

#### 问题

给定 $f, g \in \mathbb{Z}[x_1, \ldots, x_n]$，计算 $\gcd(f, g)$。多元 GCD 的挑战在于：(1) 多元多项式没有规范的"带余除法"；(2) 系数结构更复杂。

#### Brown 求值—插值策略（二元情形）

设 $f, g \in \mathbb{Z}[x, y]$，视为主变量 $x$ 的多项式，系数在 $\mathbb{Z}[y]$ 中。

**步骤**：

1. **求值**：选 $\alpha \in \mathbb{Z}$，计算 $f(x, \alpha)$ 和 $g(x, \alpha)$（一元多项式）
2. **一元 GCD**：$d_\alpha(x) = \gcd(f(x, \alpha), g(x, \alpha))$
3. **验证**：检查 $d_\alpha(x)$ 的次数是否等于真 GCD 的次数（对多个 $\alpha$ 取最小一致次数）
4. **Lagrange 插值**：对足够多的 $\alpha_1, \ldots, \alpha_k$ 计算 $d_{\alpha_i}(x)$，将每个 $x$ 系数视为 $y$ 的多项式，通过 Lagrange 插值恢复 $d(x, y)$

**关键条件**：$\alpha$ 必须是"幸运的"——即 $\gcd(f(x, \alpha), g(x, \alpha))$ 的次数等于真 GCD 的次数。若 $\alpha$ 不幸（例如使首项系数消失），需重试。

#### 多元情形的递归

$n$ 变量的 GCD 递归地减少一个变量：

1. 选 $x_n$ 的求值点 $\alpha$，计算 $f|_{x_n = \alpha}$ 和 $g|_{x_n = \alpha}$
2. 递归计算 $(n-1)$ 变量的 GCD
3. 插值恢复 $x_n$ 的依赖

#### $\mathbb{Z}$ 到 $\mathbb{Q}$ 的嵌入

在 $\mathbb{Z}[x_1, \ldots, x_n]$ 上计算 GCD 时，内容的处理至关重要：

1. 先提取关于主变量 $x_0$ 的**内容**（content in $x_0$），即所有 $x_0^k$ 系数的 GCD
2. 对本原部分进行多元 GCD
3. 最终 GCD = $\gcd(\text{cont}_x(f), \text{cont}_x(g)) \cdot \gcd(\text{pp}_x(f), \text{pp}_x(g))$

### Hensel 提升

#### 核心思想

Hensel 提升是将模 $p$ 的因式分解"提升"到模 $p^k$（或 $\mathbb{Z}$）的关键技术。设 $f \in \mathbb{Z}[x]$ 是首一多项式，且

$$f \equiv g_0 \cdot h_0 \pmod{p}, \quad \gcd(g_0, h_0) = 1 \pmod{p}$$

目标是找到 $g, h \in \mathbb{Z}[x]$ 使得 $f = g \cdot h$，且 $g \equiv g_0 \pmod{p}$，$h \equiv h_0 \pmod{p}$。

#### 线性 Hensel 提升（两因子情形）

设 $s, t$ 是 $\mathbb{F}_p[x]$ 上的 Bézout 系数：$s \cdot g_0 + t \cdot h_0 = 1 \pmod{p}$。

**迭代步骤**（从 $m = p$ 到 $m > B$）：

```
1. e ← f - g·h                    // 误差
2. if e = 0: return (g, h)        // 精确分解
3. ē ← (e / m) mod p              // 误差的"导数"
4. Δg ← (t · ē) mod g₀           // g 的修正（deg < deg g₀）
5. Δh ← (ē - Δg · h₀) / g₀      // h 的修正（精确除法）
6. g ← g + Δg · m
7. h ← h + Δh · m
8. m ← m · p
```

每步迭代将精度从 $\bmod m$ 提升到 $\bmod m \cdot p$。

#### Mignotte 界

提升到何时停止？我们需要一个因子系数大小的上界。**Landau–Mignotte 界**给出：

$$\|g\|_\infty \le 2^n \|f\|_2$$

其中 $n = \deg f$，$\|f\|_2 = \sqrt{\sum a_i^2}$ 是 2-范数。因此只需提升到 $p^k > 2 \cdot 2^n \|f\|_2$。

#### 多因子 Hensel 提升

对于 $r > 2$ 个因子，采用**逐因子剥离**策略：

```
输入：f, [g₁, g₂, …, gᵣ] mod p
输出：[G₁, G₂, …, Gᵣ] 使得 f = ∏ Gᵢ

1. 对 i = 1, …, r-1:
2.     h₀ ← g_{i+1} · g_{i+2} · … · gᵣ  (剩余因子之积)
3.     (Gᵢ, H) ← hensel_lift_pair(f_current, gᵢ, h₀, p, bound)
4.     f_current ← H
5. Gᵣ ← f_current
```

每步都是两因子的 Hensel 提升。

#### Zassenhaus 因子重组

Hensel 提升给出 $r$ 个模 $p^k$ 的因子 $\tilde{g}_1, \ldots, \tilde{g}_r$。$\mathbb{Z}[x]$ 上的真因子是这些 $\tilde{g}_i$ 的某个子集的乘积（模 $p^k$ 约化到对称范围 $(-p^k/2, p^k/2]$）。

**重组策略**：

```
1. remaining ← [g̃₁, …, g̃ᵣ]
2. result ← []
3. for size = 1, 2, …, r:
4.     for each subset S ⊂ remaining of size |S|:
5.         candidate ← primitive_part(lc(rest) · ∏_{i∈S} g̃ᵢ)
6.         reduce candidate to symmetric range mod p^k
7.         if candidate divides f:
8.             result.append(candidate)
9.             remaining ← remaining \ S
10.            rest ← f / (∏ result)
11.            break  // 从 size=1 重新开始
12. return result
```

**关键优化**：乘以当前余因子的首项系数后再取本原部分，这是 Zassenhaus 的经典技巧——真因子 $h$ 满足 $\text{lc}(\text{rest}) \cdot \prod S = c \cdot h$。

**复杂度**：最坏情况 $2^r$ 个子集（$r$ 为 mod-$p$ 因子数），但实践中因子数通常很少。

在 oCAS 中，Hensel 提升和 Zassenhaus 重组实现在 `ocas-poly/src/factor/hensel.rs`：

```rust
// 完整的 ℤ[x] 因式分解管线：
// 1. 无平方分解
// 2. 对每个无平方分量：
//    a. 选素数 p（不整除首项系数，且 f mod p 无平方）
//    b. 在 𝔽_p 上因式分解（Cantor-Zassenhaus）
//    c. Hensel 提升到 ℤ
//    d. Zassenhaus 重组
// 入口函数：DenseUnivariatePolynomial::factor()
```

### Berlekamp 算法

#### 适用场景

在小素数 $p$ 的有限域 $\mathbb{F}_p$ 上分解无平方多项式。当 $p$ 较小时（$p \le 1000$），矩阵方法比 Cantor–Zassenhaus 的随机算法更高效。

#### 理论基础

设 $f \in \mathbb{F}_p[x]$ 是无平方的 $n$ 次多项式。定义 **Frobenius 矩阵** $Q \in \mathbb{F}_p^{n \times n}$：

$$Q_{ij} = [x^j] \text{ 的系数在 } x^{ip} \bmod f \text{ 中}$$

即 $Q$ 的第 $i$ 行是 $x^{ip} \bmod f$ 的系数向量。

**关键定理**：$v \in \mathbb{F}_p^n$ 满足 $Q^T v = v$（即 $v$ 是 $Q^T - I$ 的零空间向量）当且仅当对应的多项式 $v(x) = \sum v_i x^i$ 满足

$$v(x)^p \equiv v(x) \pmod{f}$$

这意味着 $v(x)$ 在 $\mathbb{F}_p$ 中取值，且

$$f \mid \prod_{a \in \mathbb{F}_p} (v(x) - a)$$

因此 $\gcd(f, v(x) - a)$ 对某个 $a \in \mathbb{F}_p$ 给出 $f$ 的非平凡因子。

#### 算法步骤

```
输入：monic square-free f ∈ 𝔽_p[x]，deg f = n
输出：f 的不可约因子列表

1. 构造 Frobenius 矩阵 Q
2. 计算 Q^T - I 的零空间基 {v₁, v₂, …, vᵣ}
3. factors ← [f]
4. for each vⱼ (r ≥ 1):
5.     new_factors ← []
6.     for each factor g in factors:
7.         if deg(g) ≤ 1:
8.             new_factors.append(g); continue
9.         for a = 0, 1, …, p-1:
10.            d ← gcd(g, vⱼ - a)
11.            if 0 < deg(d) < deg(g):
12.                new_factors.append(d)
13.                g ← g / d
14.        new_factors.append(g)
15.    factors ← new_factors
16. return factors
```

**零空间维数**：若 $f$ 有 $r$ 个不可约因子，则零空间维数为 $r$。因此只需 $r - 1$ 个非平凡零空间向量即可完全分裂。

**复杂度**：$O(n^3)$（矩阵消元）+ $O(r \cdot n^2 \cdot p)$（GCD 分裂）。当 $p$ 很大时，$a$-循环的代价过高，应改用 Cantor–Zassenhaus。

在 oCAS 中，`berlekamp()` 实现在 `ocas-poly/src/factor/finite_field.rs`。当 $p \le 1000$ 时自动使用。

### Cantor–Zassenhaus 算法

#### 适用场景

大素数 $p$ 的有限域 $\mathbb{F}_p$ 上的无平方多项式分解。分为两个阶段：**不同次分解**（DDF）和**等次分解**（EDF）。

#### 不同次分解（DDF）

**目标**：将 $f$ 分解为 $f = g_1 \cdot g_2 \cdots g_s$，其中 $g_d$ 的每个不可约因子的次数恰好为 $d$。

**理论基础**：$\mathbb{F}_{p^d}$ 的元素恰好是 $x^{p^d} - x = 0$ 的根。因此 $f$ 的次数为 $d$ 的不可约因子恰好整除

$$\gcd(f,\, x^{p^d} - x)$$

但不整除之前的 $x^{p^{d'}} - x$（$d' < d$）。

**算法**：

```
输入：monic square-free f ∈ 𝔽_p[x]
输出：[(g₁, 1), (g₂, 2), …] 其中 gₖ 的因子次数为 k

1. current ← f,  h ← x,  degree ← 1
2. while deg(current) ≥ 2·degree:
3.     h ← h^p mod current           // Frobenius 迭代：h = x^(p^degree)
4.     g ← gcd(current, h - x)
5.     if deg(g) > 0:
6.         输出 (monic(g), degree)
7.         current ← current / g
8.         h ← h mod current
9.     degree ← degree + 1
10. if deg(current) > 0:
11.     输出 (monic(current), deg(current))
```

**Frobenius 迭代的优化**：不需要直接计算 $x^{p^d}$（指数巨大），而是逐步迭代 $h \leftarrow h^p \bmod f$，每次使用快速幂取模。

**复杂度**：$O(n^2 \log p)$（每步 Frobenius 迭代），总共 $O(n^3 \log p)$。

#### 等次分解（EDF）

**目标**：将 $f = g_1 g_2 \cdots g_r$（每个 $g_i$ 不可约且次数为 $d$）完全分裂。

**奇特征情形**（$p > 2$）：

利用 $\mathbb{F}_{p^d}$ 中恰好一半元素是平方剩余的事实。对随机多项式 $a$，计算

$$b = a^{(p^d - 1)/2} \bmod f$$

则 $b$ 在每个不可约因子 $g_i$ 上的取值要么是 $1$（$a$ 是该因子域中的平方剩余），要么是 $-1$（非剩余）。因此

$$\gcd(f,\, b - 1) = \prod_{i:\, a \text{ 是 } g_i \text{ 中的平方剩余}} g_i$$

给出 $f$ 的非平凡分裂（除非所有 $g_i$ 的取值一致，概率为 $2^{1-r}$）。

**特征 2 情形**（$p = 2$）：

$b - 1$ 的技巧不适用（$1 = -1$）。改用**迹映射**（trace map）：

$$T(a) = a + a^2 + a^{2^2} + \cdots + a^{2^{d-1}} \bmod f$$

$T(a)$ 在每个 $\mathbb{F}_{2^d}$ 上取值为 $\mathbb{F}_2$ 中的元素（迹的性质），且对随机 $a$ 的取值在不同因子上独立。因此 $\gcd(f, T(a))$ 给出非平凡分裂。

```
输入：f（DDF 输出的 d-次因子之积）
输出：f 的不可约因子列表

1. factors ← [f]
2. while 存在 deg > d 的因子:
3.     for each factor g with deg(g) > d:
4.         选取随机多项式 a（deg < deg(g)）
5.         if p = 2:
6.             b ← T(a) = Σᵢ₌₀^{d-1} a^{2^i} mod g
7.         else:
8.             b ← a^{(p^d-1)/2} mod g
9.         d₁ ← gcd(g, b - 1)   // 或 gcd(g, b) for char 2
10.        if 0 < deg(d₁) < deg(g):
11.            replace g by d₁ and g/d₁ in factors
12. return factors
```

在 oCAS 中，`cantor_zassenhaus()` 实现在 `ocas-poly/src/factor/finite_field.rs`，顶层入口 `factor_over_finite_field()` 根据 $p$ 的大小自动选择 Berlekamp（$p \le 1000$）或 Cantor–Zassenhaus。

### Wang EEZ 多元因式分解

#### 问题

给定 $f \in \mathbb{Z}[x_1, \ldots, x_n]$（或 $\mathbb{F}_p[x_1, \ldots, x_n]$），分解为不可约因子。

#### 策略概述

Wang（1978）的 EEZ（Evaluation, Exact division, Zassenhaus）算法将一元 Hensel 提升推广到多元：

1. **求值**：在辅助变量上选求值点，得到一元像 $f(x_1, \alpha_2, \ldots, \alpha_n)$
2. **一元分解**：在 $\mathbb{Z}$（或 $\mathbb{F}_p$）上分解一元像
3. **逐变量 Hensel 提升**：一个变量一个变量地提升回多元分解
4. **Zassenhaus 重组**：组合提升后的因子

#### Wang 首项系数预处理

多元因式分解的难点之一是非常数首项系数。设 $f$ 的主变量 $x_0$ 的首项系数为 $\ell(x_1, \ldots, x_n)$。若 $\ell$ 不是常数，每个因子 $f_i$ 的首项系数 $\ell_i$ 也不一定是常数，且 $\prod \ell_i = \ell$。

**Wang 的贪心分配**：

1. 将 $\ell$ 分解为不可约因子的幂之积：$\ell = \prod g_j^{e_j}$
2. 在求值点 $\alpha = (\alpha_2, \ldots, \alpha_n)$ 处，$\ell$ 的像分解为一元因子的像
3. 将 $\ell$ 的非平凡因子 $g_j$ 分配给对应的 $f_i$，使得 $g_j(\alpha) = \text{lc}(u_i)$（$u_i$ 是一元像因子）
4. 通过逐对互素条件（$\alpha_j = |g_j(\alpha)| > 1$）确保分配的一致性

#### 逐变量 Hensel 提升

设 $f$ 有 $n$ 个变量，一元像因子为 $u_1, \ldots, u_r$。

**提升变量 $x_k$**（$k = 1, 2, \ldots, n-1$，按变量序号逐个提升）：

1. 对当前因子 $F_i^{(k-1)}$（在变量 $x_1, \ldots, x_{k-1}$ 中），求值 $x_k = \alpha_k$ 得到 $u_i$
2. 计算 Bézout 系数 $b_i$ 使得 $\sum b_i \prod_{j \ne i} u_j = 1$
3. 对 $t = 1, 2, \ldots$ 逐次求解**多元 Diophantine 方程**：

$$\sum_i \sigma_i \cdot \prod_{j \ne i} F_j = e_t$$

其中 $e_t$ 是当前误差 $f - \prod F_i$ 的 Taylor 展开第 $t$ 项

4. 修正 $F_i \leftarrow F_i + \sigma_i \cdot (x_k - \alpha_k)^t$

#### 多元 Diophantine 方程求解器

提升变量 $x_k$ 时需要求解的方程形式为：

$$\sum_{i=1}^{r} \sigma_i \cdot \prod_{j \ne i} u_j = e \pmod{(x_k - \alpha_k)^t}$$

其中 $\deg_{x_0}(\sigma_i) < \deg_{x_0}(u_i)$。这是一个线性方程组（关于 $\sigma_i$ 的系数），可以通过以下方式递归求解：

- 当 $k = 1$（单变量）：扩展 Euclid 算法
- 当 $k > 1$：递归到 $k-1$ 个变量，在 $\alpha_k$ 处求值，求解更小的系统，再插值

#### p-adic 系数 Hensel 提升

对于 $\mathbb{Z}[x_1, \ldots, x_n]$ 上的非常数首项系数情形，Wang EEZ 之后还需要**p-adic 系数提升**：

1. 将多元因子的系数模 $p$ 得到 $\mathbb{F}_p$ 上的骨架
2. 迭代求解 mod-$p$ Diophantine 方程，逐步将系数从 $\bmod p$ 提升到 $\bmod p^k$
3. 直到误差为零或 $p^k$ 超过系数界（Gelfond 界）

**Gelfond 系数界**：

$$B = \left(\sqrt{\prod_v (d_v + 1) \cdot 2^{2 \sum d_v - n}} + 1\right) \cdot \|f\|_\infty \cdot |\text{lc}(f)|$$

其中 $d_v$ 是变量 $v$ 的次数，$n$ 是变量数。

#### 稀疏 Diophantine 求解

当 Diophantine 修正项具有稀疏结构时，可以避免密集递归求解。oCAS 实现了**骨架插值**（skeleton interpolation）策略：

1. 从误差 $e$ 的项中提取骨架（skeleton）——修正项的可能指数模式
2. 在多个随机基点上求值，得到单变量 Diophantine 方程
3. 通过 Vandermonde 系统插值系数
4. 验证插值结果

这在因子具有许多变量但每个变量的度很低时特别有效。

在 oCAS 中，通用 EEZ 算法实现在 `ocas-poly/src/factor/eez.rs`，二元特化版本在 `ocas-poly/src/factor/multivariate.rs`（目前要求主变量首项系数为常数）：

```rust
// ℤ 上的 n 元因式分解入口（eez.rs）：
pub fn multivariate_factor_z(f: &ZmPoly) -> Vec<(ZmPoly, usize)>

// 𝔽_p 上的 n 元因式分解入口（eez.rs）：
pub fn multivariate_factor_fp(f: &FpMPoly) -> Vec<(FpMPoly, usize)>

// 二元特化入口（multivariate.rs）：
pub fn bivariate_factor_z(f: &ZMPoly, x_var: usize, y_var: usize) -> Vec<(ZMPoly, usize)>
pub fn bivariate_factor_fp(f: &FpMPoly, x_var: usize, y_var: usize) -> Vec<(FpMPoly, usize)>
```

### Trager 代数数域因式分解

#### 问题

给定 $f \in K[x]$，其中 $K = \mathbb{Q}(\alpha)$ 是代数数域（$\alpha$ 是极小多项式 $m(\alpha) = 0$ 的根），分解 $f$ 为 $K[x]$ 中的不可约因子。

#### 范数下降

Trager 算法的核心思想是通过**范数**将问题从 $K[x]$ 降到 $\mathbb{Q}[x]$：

$$N(f) = \text{Res}_\alpha(m(\alpha),\, f(x, \alpha))$$

范数 $N(f) \in \mathbb{Q}[x]$ 的次数为 $\deg_x(f) \cdot [K:\mathbb{Q}]$，且若 $f$ 在 $K[x]$ 中可约，则 $N(f)$ 在 $\mathbb{Q}[x]$ 中可约。

#### 求值—插值计算范数

直接构造 Sylvester 矩阵计算结式代价高昂。oCAS 使用**求值—插值**：

1. 在 $\deg_x(f) \cdot [K:\mathbb{Q}] + 1$ 个有理点 $x_j$ 处计算标量结式 $\text{Res}_\alpha(m, f(x_j, \alpha))$
2. 每个标量结式是 $\mathbb{Q}(\alpha)$ 中元素的范数，通过在 $[K:\mathbb{Q}]$ 个点处求值 $m$ 并取乘积得到
3. 通过 Newton 差商插值恢复 $N(f) \in \mathbb{Q}[x]$

#### Trager 平移

若 $N(f)$ 有重因子，$f$ 的因子信息会丢失。**Trager 平移**通过替换 $x \mapsto x - s\alpha$（$s \ge 0$）使范数无平方：

$$N(f(x - s\alpha)) \text{ 无平方}$$

这样的 $s$ 总是存在的（坏的 $s$ 只有有限多个），实践中 $s = 0$ 或很小的值即可。

#### 分解步骤

```
输入：f ∈ K[x]（无平方）
输出：f 在 K[x] 中的不可约因子

1. for s = 0, 1, 2, …:
2.     g ← f(x - s·α)
3.     N ← Res_α(m, g(x, α))  via 求值-插值
4.     if N 无平方: break
5. 
6. N₁, …, Nₖ ← factor_over_Q(N)   // ℚ 上因式分解
7. for each Nᵢ:
8.     N̂ᵢ ← embed_Q_to_K(Nᵢ)
9.     hᵢ ← gcd_K(g, N̂ᵢ)          // K 上的 GCD
10.    output ← compose_linear(hᵢ, +s·α)   // 逆平移
11. return [monic(output₁), …]
```

#### 代数数域上的模 GCD

步骤 9 中的 $K[x]$ GCD 通过模方法高效计算（`gcd_anf`）：

1. 选素数 $p$ 使得 $m$ 在 $\mathbb{F}_p$ 上不可约（确保 $\mathbb{F}_p[\alpha]/(m) \cong \text{GF}(p^d)$）
2. 将 $a, b \in K[x]$ 映射到 $\text{GF}(p^d)[x]$
3. 计算 $\text{GF}(p^d)[x]$ 上的 monic GCD
4. CRT 合并 + 有理重构 + 试除验证

不幸素数（模 GCD 次数偏大）被丢弃，最多尝试 1000 个素数后回退到密集 Euclid GCD。

在 oCAS 中，Trager 算法实现在 `ocas-poly/src/factor/algebraic.rs`：

```rust
// K = ℚ(α) 上的因式分解入口：
impl DenseUnivariatePolynomial<AlgebraicNumberField> {
    pub fn factor(&self) -> Factors<AlgebraicNumberField>
}

// 示例：x² - 2 在 ℚ(√2) 上分解为 (x - √2)(x + √2)
let field = AlgebraicNumberField::new(
    RationalDomain,
    vec![Rational::new(-2, 1), Rational::new(0, 1), Rational::new(1, 1)],
);
let f = DenseUnivariatePolynomial::from_coeffs(field.clone(), vec![
    field.from_base(Rational::new(-2, 1)),
    field.zero(),
    field.one(),
]);
let factors = f.factor();
assert_eq!(factors.len(), 2);  // (x - √2) 和 (x + √2)
```

## 在 oCAS 中的实现

### 文件映射

| 源文件 | 算法 | 说明 |
|---|---|---|
| `ocas-poly/src/dense.rs` | Euclidean GCD、Karatsuba 乘法 | 一元多项式的核心数据结构 |
| `ocas-poly/src/gcd/modular.rs` | Brown 模 GCD | $\mathbb{Z}[x]$ 上的高效 GCD |
| `ocas-poly/src/factor/mod.rs` | 顶层入口、Yun 无平方分解 | `factor()` 和 `square_free_factorization()` |
| `ocas-poly/src/factor/finite_field.rs` | Berlekamp、Cantor–Zassenhaus（DDF + EDF） | $\mathbb{F}_p[x]$ 因式分解 |
| `ocas-poly/src/factor/hensel.rs` | Hensel 提升、Zassenhaus 重组、Mignotte 界 | $\mathbb{Z}[x]$ 因式分解管线 |
| `ocas-poly/src/factor/multivariate.rs` | 二元因式分解（Wang EEZ 特化；主变量首项系数需为常数） | $\mathbb{Z}[x,y]$ 和 $\mathbb{F}_p[x,y]$ |
| `ocas-poly/src/factor/eez.rs` | 通用多元 EEZ Hensel 提升 | $\mathbb{Z}[x_1,\ldots,x_n]$ 和 $\mathbb{F}_p[x_1,\ldots,x_n]$ |
| `ocas-poly/src/factor/algebraic.rs` | Trager 算法（代数数域因式分解） | $\mathbb{Q}(\alpha)[x]$ |
| `ocas-poly/src/multivariate_gcd.rs` | 多元 GCD（求值—插值） | 辅助多元因式分解 |

### 算法选择策略

oCAS 根据输入自动选择最优算法：

```
factor(f ∈ ℤ[x]):
  1. 无平方分解（Yun 算法）
  2. 对每个无平方分量 g:
     a. if deg(g) ≤ 1: 直接返回
     b. 选素数 p（不整除 lc(g)，g mod p 无平方）
     c. 在 𝔽_p 上分解（Cantor–Zassenhaus；通用的 factor_over_finite_field 入口则在 p ≤ 1000 时自动改用 Berlekamp）
     d. if mod-p 因子数 = 1: g 不可约，返回
     e. Hensel 提升（Mignotte 界决定精度）
     f. Zassenhaus 因子重组
     g. 非首一情况：首项系数变换 a^{d-1}·f(x/a)

factor(f ∈ 𝔽_p[x]):
  1. 提取首项系数
  2. 无平方分解（Musser/Bernardin，处理特征 p）
  3. 对每个无平方分量：Berlekamp（小 p）或 Cantor-Zassenhaus

factor(f ∈ ℚ(α)[x]):
  1. Yun 无平方（使用模 GCD gcd_anf）
  2. 对每个分量：Trager 范数 + ℚ 分解 + K-GCD 恢复

factor(f ∈ ℤ[x₁,…,xₙ]):
  1. 提取内容（主变量 x₀）
  2. 无平方分解
  3. Wang 首项系数预处理
  4. 选求值点，一元分解
  5. EEZ 逐变量 Hensel 提升
  6. p-adic 系数提升（非常数首项系数时）
  7. Zassenhaus 重组
```

### 性能特征

| 算法 | 时间复杂度 | 空间复杂度 | 适用场景 |
|---|---|---|---|
| Euclidean GCD（PRS） | $O(n^2)$ | $O(n)$ | 低度（$\lesssim 16$） |
| Brown 模 GCD | $O(n^2 \cdot k)$ | $O(n \cdot k)$ | 高度 $\mathbb{Z}[x]$，$k$ = 素数个数 |
| Hensel 提升 | $O(n^2 \cdot \log B)$ | $O(n \cdot r)$ | $r$ 因子，$B$ = Mignotte 界 |
| Zassenhaus 重组 | $O(2^r \cdot n)$ | $O(n \cdot r)$ | $r$ 个 mod-$p$ 因子 |
| Berlekamp | $O(n^3 + r \cdot n^2 p)$ | $O(n^2)$ | 小 $p$ |
| Cantor–Zassenhaus | $O(n^3 \log p)$ | $O(n)$ | 大 $p$ |
| Wang EEZ | 指数（变量数） | 多元多项式 | 多元分解 |
| Trager | $O(n^3 [K:\mathbb{Q}]^2)$ | $O(n \cdot [K:\mathbb{Q}])$ | 代数数域 |

## 参考文献

- **[Gathen–Gerhard]** J. von zur Gathen and J. Gerhard, *Modern Computer Algebra*, 3rd ed., Cambridge University Press, 2013. Chapters 14–15 (factoring), 18 (GCD).
- **[Brown]** W. S. Brown, "On Euclid's Algorithm and the Computation of Polynomial Greatest Common Divisors," *J. ACM*, 18(4):478–504, 1971.
- **[Zassenhaus]** H. Zassenhaus, "On Hensel Factorization I," *J. Number Theory*, 1(3):291–311, 1969.
- **[Wang]** P. S. Wang, "An Improved Multivariate Polynomial Factoring Algorithm," *Math. Comp.*, 32(144):1215–1231, 1978.
- **[Berlekamp]** E. R. Berlekamp, "Factoring Polynomials over Large Finite Fields," *Math. Comp.*, 24(111):713–735, 1970.
- **[Cantor–Zassenhaus]** D. G. Cantor and H. Zassenhaus, "A New Algorithm for Factoring Polynomials over Finite Fields," *Math. Comp.*, 36(154):587–592, 1981.
- **[Trager]** B. M. Trager, "Algebraic Factoring and Rational Function Integration," *Proc. SYMSAC '76*, pp. 219–226, 1976.
- **[Mignotte]** M. Mignotte, *Mathematics for Computer Algebra*, Springer, 1992.
- **[Geddes–Czapor–Labahn]** K. O. Geddes, S. R. Czapor, and G. Labahn, *Algorithms for Computer Algebra*, Kluwer, 1992. Chapters 6 (multivariate factoring), 7 (GCD).
- **[Knuth]** D. E. Knuth, *The Art of Computer Programming*, Vol. 2, §4.6.2, Addison-Wesley, 1997.
