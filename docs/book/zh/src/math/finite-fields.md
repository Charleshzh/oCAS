# 基础：有限域与模算术

## 前提知识

- 集合论与基本代数结构（群、环、域的定义）
- 整数除法与余数
- 多项式的基本运算

建议先阅读：[多项式代数](./polynomial-algebra.md)、[线性代数](./linear-algebra.md)。

---

## 基础概念

### 模算术

给定正整数 $n \geq 2$，整数 $a$ 与 $b$ **模 $n$ 同余**（记为 $a \equiv b \pmod{n}$）当且仅当 $n \mid (a - b)$。模运算将整数划分为 $n$ 个等价类：

$$
\mathbb{Z}/n\mathbb{Z} = \{\overline{0},\, \overline{1},\, \dots,\, \overline{n-1}\}
$$

其中 $\overline{a} = \{a + kn \mid k \in \mathbb{Z}\}$。加法和乘法按代表元运算后再取模：

$$
\overline{a} + \overline{b} = \overline{a + b}, \quad \overline{a} \cdot \overline{b} = \overline{a \cdot b}
$$

这使得 $\mathbb{Z}/n\mathbb{Z}$ 成为一个交换环。

### 模逆与扩展 Euclid 算法

元素 $\overline{a} \in \mathbb{Z}/n\mathbb{Z}$ 有乘法逆元当且仅当 $\gcd(a, n) = 1$。此时存在整数 $x, y$ 满足 Bézout 等式：

$$
ax + ny = \gcd(a, n) = 1
$$

取模 $n$ 即得 $ax \equiv 1 \pmod{n}$，故 $\overline{x} = \overline{a}^{-1}$。

**扩展 Euclid 算法**通过迭代求解此等式：

1. 初始化：$(r_0, r_1) = (a, n)$，$(s_0, s_1) = (1, 0)$，$(t_0, t_1) = (0, 1)$。
2. 每步计算 $q = \lfloor r_{i-1}/r_i \rfloor$，更新：

$$
\begin{pmatrix} r_{i+1} \\ s_{i+1} \\ t_{i+1} \end{pmatrix} = \begin{pmatrix} r_{i-1} - q \cdot r_i \\ s_{i-1} - q \cdot s_i \\ t_{i-1} - q \cdot t_i \end{pmatrix}
$$

3. 当 $r_{i+1} = 0$ 时终止，此时 $r_i = \gcd(a, n)$，$s_i \cdot a + t_i \cdot n = r_i$。

### 素域 $\mathbb{Z}/p\mathbb{Z}$

当 $n = p$ 为**素数**时，$\mathbb{Z}/p\mathbb{Z}$ 的结构发生根本变化：每个非零元素都与 $p$ 互素，因此都有逆元。这使得 $\mathbb{Z}/p\mathbb{Z}$ 成为一个**域**（field），记为 $\mathbb{F}_p$。

**定理（Fermat 小定理）**：设 $p$ 为素数，$a \not\equiv 0 \pmod{p}$，则

$$
a^{p-1} \equiv 1 \pmod{p}
$$

由此立即得到 $a^{p-2} \equiv a^{-1} \pmod{p}$，这为有限域求逆提供了一种**无需扩展 Euclid 算法**的实现方式——只需一次模幂运算。当使用 GMP 后端的 `pow_mod` 时，这比逐次除法快得多。

**域的基本性质**：

- 加法群 $(\mathbb{F}_p, +)$ 是 $p$ 阶循环群，同构于 $(\mathbb{Z}/p\mathbb{Z}, +)$。
- 乘法群 $(\mathbb{F}_p^*, \times)$ 是 $p-1$ 阶循环群——存在**原根**（primitive root）$g$ 使得 $\{g^0, g^1, \dots, g^{p-2}\} = \mathbb{F}_p^*$。
- **特征**（characteristic）为 $p$：$p \cdot 1 = \underbrace{1 + \cdots + 1}_{p} = 0$，且 $p$ 是使此式成立的最小正整数。

---

## 核心理论

### 有限域 $\mathbb{F}_{p^d}$ 的构造

素域 $\mathbb{F}_p$ 并非唯一的有限域。**有限域分类定理**指出：对每个素数幂 $q = p^d$（$d \geq 1$），存在唯一的（在同构意义下）$q$ 元有限域 $\mathbb{F}_q$，且所有有限域都是这种形式。

构造 $\mathbb{F}_{p^d}$（$d > 1$）的标准方法是**不可约多项式扩张**：

1. 选取 $\mathbb{F}_p[x]$ 中一个 $d$ 次**不可约多项式** $m(x)$。
2. 定义 $\mathbb{F}_{p^d} = \mathbb{F}_p[x] / (m(x))$——多项式环对理想 $(m(x))$ 的商环。
3. 元素是次数 $< d$ 的多项式 $a_0 + a_1 \alpha + \cdots + a_{d-1}\alpha^{d-1}$，其中 $\alpha = \overline{x}$ 是 $m$ 的根。
4. 加法按系数逐项模 $p$；乘法先做多项式乘法，再对 $m(x)$ 取模。

**不可约多项式的存在性**：$\mathbb{F}_p[x]$ 中**首一**（monic）$d$ 次不可约多项式的个数为

$$
N_p(d) = \frac{1}{d}\sum_{k \mid d} \mu(k) \cdot p^{d/k}
$$

其中 $\mu$ 是 Möbius 函数。这保证了对任意 $d \geq 1$，不可约多项式都存在。

**例**：$\mathbb{F}_4 = \mathbb{F}_2[x]/(x^2 + x + 1)$。设 $\alpha$ 为 $x^2 + x + 1$ 的根，则 $\alpha^2 = \alpha + 1$（在 $\mathbb{F}_2$ 中 $-1 = 1$）。四个元素为 $\{0, 1, \alpha, \alpha+1\}$。

### 循环乘法群

**定理**：有限域的乘法群 $\mathbb{F}_{p^d}^*$ 是循环群，阶为 $p^d - 1$。

这意味着存在**本原元**（primitive element）$g \in \mathbb{F}_{p^d}^*$ 使得

$$
\mathbb{F}_{p^d}^* = \{g^0, g^1, \dots, g^{p^d - 2}\}
$$

循环性在密码学和编码理论中有深远应用：离散对数问题（DLP）、伪随机数生成、纠错码的构造都依赖于此。

### 特征与素子域

有限域 $\mathbb{F}_{p^d}$ 的**特征**为 $p$。由特征定义：

- $p \cdot a = 0$ 对所有 $a \in \mathbb{F}_{p^d}$ 成立。
- $\mathbb{F}_p$ 是 $\mathbb{F}_{p^d}$ 的**素子域**（prime subfield），即包含 $1$ 的最小子域。
- $\mathbb{F}_{p^d}$ 作为 $\mathbb{F}_p$-向量空间，维数为 $d$。

Frobenius 自同态 $\varphi: x \mapsto x^p$ 是 $\mathbb{F}_{p^d}$ 上的重要结构映射。其 $d$ 次迭代 $\varphi^d = \mathrm{id}$，不动点恰好是 $\mathbb{F}_p$ 中的元素。

---

## 在 oCAS 中的实现

oCAS 的有限域实现在 `ocas-domain/src/finite_field.rs` 中，支持素域 $\mathbb{Z}/p\mathbb{Z}$。扩展域 $\mathbb{F}_{p^d}$ 通过 `AlgebraicExtension<FiniteField>` 实现（见 `ocas-domain/src/algebraic.rs`）。

### 核心结构体

**`FiniteField`** 表示一个素域 $\mathbb{F}_p$：

```rust
pub struct FiniteField {
    prime: BigInt,
    prime_minus_two: BigInt,  // 缓存 p-2，用于 Fermat 求逆
    #[cfg(feature = "gmp")]
    prime_gmp: rug::Integer,  // GMP 后端加速模运算
    #[cfg(feature = "gmp")]
    prime_minus_two_gmp: rug::Integer,
}
```

- `prime` 字段存储模数 $p$，使用 `num-bigint` 的任意精度整数。
- `prime_minus_two` 在构造时预计算，避免每次求逆时重复计算。
- 启用 `gmp` feature 时，额外缓存 GMP 表示以利用 `rug` 的原生模幂。

**`FiniteFieldElement`** 表示域中元素：

```rust
pub struct FiniteFieldElement {
    value: BigInt,  // 始终保持在 [0, p-1] 范围内
}
```

元素始终以**规范代表元**（canonical representative）存储——即 $[0, p-1]$ 范围内的非负整数。这保证了语义相等即结构相等（`PartialEq` 直接比较值），无需每次比较时取模。

### 构造与元素创建

```rust
// 创建 F_7
let f = FiniteField::new(BigInt::from(7));

// 从任意整数创建元素（自动归约到 [0, p-1]）
let a = f.element(10);   // a = 10 mod 7 = 3
let b = f.element(-3);   // b = -3 mod 7 = 4
```

`element()` 方法使用 `mod_floor`（而非简单的 `%`），确保负数也被正确映射到 $[0, p-1]$。

### 算术运算

所有运算通过 `Domain` trait 统一接口：

| 运算 | 方法 | 实现策略 |
|---|---|---|
| 加法 | `f.add(&a, &b)` | $(a + b) \bmod p$ |
| 减法 | `f.sub(&a, &b)` | $(a - b) \bmod p$（使用 `mod_floor` 处理负数） |
| 乘法 | `f.mul(&a, &b)` | $(a \cdot b) \bmod p$ |
| 取反 | `f.neg(&a)` | $(-a) \bmod p$ |
| 求逆 | `f.inv(&a)` | $a^{p-2} \bmod p$（Fermat 小定理）；$a = 0$ 时返回 `None` |
| 除法 | `f.div(&a, &b)` | $a \cdot b^{-1}$；$b = 0$ 时返回 `None` |
| 幂运算 | `f.pow(&a, n)` | 快速模幂（二进制取幂） |

**求逆实现细节**：

```rust
fn inv(&self, a: &FiniteFieldElement) -> Option<FiniteFieldElement> {
    if a.value.is_zero() {
        return None;  // 零没有逆元
    }
    // Fermat: a^(p-2) ≡ a^{-1} (mod p)
    Some(self.normalize(a.value.modpow(&self.prime_minus_two, &self.prime)))
}
```

使用 Fermat 小定理而非扩展 Euclid 算法，因为 `modpow` 的二进制取幂只需 $O(\log p)$ 次乘法，且实现简单。当启用 GMP 后端时，`rug::Integer::pow_mod` 使用底层 GMP 的高度优化实现。

### EuclideanDomain 实现

有限域同时实现了 `EuclideanDomain` trait，但域上的欧几里得除法是退化的：

```rust
// 域上除法总是精确的，余数恒为零
fn div_rem(&self, a, b) -> Option<(Self::Element, Self::Element)> {
    self.div(a, b).map(|q| (q, self.zero()))
}

// 域上的 GCD：两者均为零时返回 0，否则返回 1
fn gcd(&self, a, b) -> Self::Element {
    if self.is_zero(a) && self.is_zero(b) {
        self.zero()
    } else {
        self.one()
    }
}
```

这是因为域中每个非零元素都是**单位**（unit），因此 GCD 的概念退化为平凡情况。

### 扩展域 $\mathbb{F}_{p^d}$

通过 `AlgebraicExtension<FiniteField>` 构造：

```rust
use ocas_domain::{AlgebraicExtension, Domain, FiniteField};
use num_bigint::BigInt;

// 构造 F_4 = F_2[x]/(x^2 + x + 1)
let f2 = FiniteField::new(BigInt::from(2));
let f4 = AlgebraicExtension::new(
    f2,
    vec![
        f2.element(1),  // 常数项
        f2.element(1),  // x 系数
        f2.element(1),  // x^2 系数（首项）
    ],
);
let alpha = f4.alpha();  // α = x mod (x^2+x+1)
// α^2 = α + 1（在 F_2 中）
let alpha_sq = f4.mul(&alpha, &alpha);
```

`AlgebraicElement<E>` 内部用系数向量表示，按升幂排列、尾部零裁剪：

```rust
pub struct AlgebraicElement<E> {
    coeffs: Vec<E>,  // [a_0, a_1, ..., a_{d-1}] 表示 a_0 + a_1*α + ...
}
```

求逆使用**扩展 Euclid 算法**（在基域的稠密多项式上自包含实现），求得 $s(x) \cdot a(x) + t(x) \cdot m(x) = 1$，则 $s(x) \bmod m(x) = a(x)^{-1}$。

---

## 进阶话题

### 数论变换（NTT）与快速多项式乘法

NTT 是 FFT 在有限域上的类比。标准 FFT 在复数域 $\mathbb{C}$ 上使用单位根 $\omega_n = e^{2\pi i/n}$，NTT 则在 $\mathbb{F}_p$ 上使用**模 $p$ 的 $n$ 次单位根**。

**存在条件**：$\mathbb{F}_p$ 中存在 $n$ 次本原单位根当且仅当 $n \mid (p - 1)$。这是因为 $\mathbb{F}_p^*$ 是 $p-1$ 阶循环群，$n$ 次单位根存在 $\iff$ $n$ 整除群的阶。

**算法流程**（给定多项式 $f, g \in \mathbb{F}_p[x]$，次数各 $< n$）：

1. **补零**：将 $f, g$ 的系数向量补到长度 $N = 2^{\lceil \log_2(2n) \rceil}$。
2. **正变换**：$\hat{f} = \mathrm{NTT}(f)$，$\hat{g} = \mathrm{NTT}(g)$。
3. **逐点乘**：$\hat{h}_i = \hat{f}_i \cdot \hat{g}_i$（逐分量模 $p$）。
4. **逆变换**：$h = \mathrm{NTT}^{-1}(\hat{h})$，最后乘以 $N^{-1} \bmod p$。

总复杂度 $O(N \log N)$，对比朴素卷积的 $O(n^2)$ 和 Karatsuba 的 $O(n^{1.585})$。

**oCAS 中的 NTT 实现**（`ocas-poly/src/ntt.rs`）：

- **算法**：Cooley-Tukey 基 2 DIT（时间抽取），bit-reversal 置换。
- **Montgomery 算术**：所有中间值以 Montgomery 形式 $x \cdot R \bmod p$（$R = 2^{64}$）存储，避免昂贵的 `u128 % p` 运算。入口和出口各做一次转换。
- **触发条件**：当多项式系数个数 $\geq 256$（`NTT_THRESHOLD`，即次数 $\geq 255$）且 $p$ 为 NTT-friendly（$p-1$ 含足够大的 2-幂因子）时自动激活，否则回退到 Karatsuba。
- **限制**：当前仅支持能放入 `u64` 的素数（$p < 2^{64}$，`try_ntt_mul_fp` 中检查）；Montgomery 约减使用 128 位中间值，$p \leq 2^{63}$ 时可保证不溢出。

```rust
// NTT-friendly 检查
pub fn is_ntt_friendly(p: u64, n: usize) -> bool {
    // n 必须整除 p-1
    if n == 0 { return true; }
    let pm1 = p - 1;
    pm1 % (n as u64) == 0
}
```

**NTT-friendly 素数的构造**：实践中常用形如 $p = k \cdot 2^m + 1$ 的素数（如 NTT Prime $998244353 = 119 \cdot 2^{23} + 1$），保证 $2^{23} \mid (p-1)$，从而可处理长度达 $2^{23}$ 的变换。

**Montgomery 乘法**的核心思想是用乘法和右移代替取模：

$$
\mathrm{MontMul}(a, b) = a \cdot b \cdot R^{-1} \bmod p
$$

其中 $R = 2^{64}$，$R^{-1} \bmod p$ 和 $p' = -p^{-1} \bmod R$ 预计算。乘法步骤为：

1. $t = a \cdot b$（128 位中间值）
2. $m = (t \bmod R) \cdot p' \bmod R$
3. $u = (t + m \cdot p) / R$
4. 若 $u \geq p$ 则 $u = u - p$

全程仅用加法、乘法和移位——无昂贵的除法。

### NTT 在多项式 GCD 中的应用

NTT 不仅加速乘法，还间接加速了多项式 GCD 和因式分解。在 Hensel 提升过程中，大量小系数多项式的乘法是瓶颈，NTT 将其从 $O(n^2)$ 降为 $O(n \log n)$。oCAS 的 `DenseUnivariatePolynomial<FiniteField>::mul_ntt` 在系数个数 $\geq 256$ 且素数 NTT-friendly 时自动启用 NTT（否则回退 Karatsuba/Schoolbook）；通用 `mul` 目前仍走 Karatsuba/Schoolbook 路径。

---

## 参考文献

1. **Shoup, V.** *A Computational Introduction to Number Theory and Algebra.* Cambridge University Press, 2nd edition, 2009. — 第 4–5 章涵盖有限域构造、循环群和不可约多项式。
2. **Lidl, R. & Niederreiter, H.** *Introduction to Finite Fields and Their Applications.* Cambridge University Press, 1994. — 有限域理论的经典教材。
3. **Gathen, J. von zur & Gerhard, J.** *Modern Computer Algebra.* Cambridge University Press, 3rd edition, 2013. — 第 8 章（快速多项式乘法与 NTT）和第 14 章（有限域算术）。
4. **Menezes, A., van Oorschot, P. & Vanstone, S.** *Handbook of Applied Cryptography.* CRC Press, 1996. — 第 2 章涵盖模算术和扩展 Euclid 算法。
