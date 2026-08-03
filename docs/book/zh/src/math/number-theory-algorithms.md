# 进阶：数论算法

本章系统讲解 oCAS 中实现的核心数论算法：素性判定、整数分解、离散对数、中国剩余定理、二次剩余以及有理重构。每个算法从数学原理出发，逐步深入到 oCAS 的具体实现细节。

## 前提知识

阅读本章前，读者应熟悉以下概念：

- **素数与整除**：素数定义、最大公因数、Euclid 算法、扩展 Euclid 算法
- **模算术**：同余、模逆（扩展 Euclid 求 $a^{-1} \bmod m$）、快速模幂（binary exponentiation）
- **Euler φ 函数**：$\varphi(n) = n \prod_{p \mid n}(1 - 1/p)$，计数 $[1, n]$ 中与 $n$ 互素的整数个数
- **积性函数**：$f(mn) = f(m)f(n)$ 当 $\gcd(m, n) = 1$；完全积性则对所有 $m, n$ 成立
- **Fermat 小定理**：$p$ 素数且 $\gcd(a, p) = 1$ 时 $a^{p-1} \equiv 1 \pmod{p}$
- **群论基础**：循环群、元素的阶、生成元（原根）

**推荐阅读**：Shoup《A Computational Introduction to Number Theory and Algebra》Ch.1–9。

## 基础概念

### 素性判定问题

给定正整数 $n$，判断 $n$ 是否为素数。这是一个判定问题，但确定性算法（如 AKS）的实际速度远不如概率性算法。现代实践采用**概率素数测试**（probable prime test）：合数通过测试的概率可压到极低，但无法绝对排除。

### 整数分解问题

给定合数 $n$，找到非平凡因子 $d$（$1 < d < n$）。重复递归直到所有因子均为素数，得到 $n = p_1^{e_1} \cdots p_k^{e_k}$。分解的难度是 RSA 等公钥密码安全性的基础假设。

### 离散对数问题

在群 $\mathbb{Z}_p^*$ 中，给定 $g$ 和 $h$，求 $x$ 使得 $g^x \equiv h \pmod{p}$。当 $p - 1$ 光滑（只有小素因子）时，Pohlig–Hellman 算法可高效求解。

### 二次剩余

对奇素数 $p$ 和整数 $a$，若存在 $x$ 满足 $x^2 \equiv a \pmod{p}$，则称 $a$ 是模 $p$ 的**二次剩余**（quadratic residue）。Legendre 符号 $\left(\frac{a}{p}\right) = 1$ 表示是，$-1$ 表示非，$0$ 表示 $p \mid a$。

### 积性数论函数

| 函数 | 定义 | 积性 |
|---|---|---|
| $\varphi(n)$ | $\#\{k \in [1,n] : \gcd(k,n)=1\}$ | 是 |
| $\mu(n)$ | $1$（$n=1$），$(-1)^k$（$n$ 恰好 $k$ 个不同素因子），$0$（$n$ 有平方因子） | 是 |
| $\tau(n)$ | 正因子个数 $= \prod(e_i + 1)$ | 是 |
| $\sigma_k(n)$ | $\sum_{d \mid n} d^k$ | 是 |
| $\lambda(n)$ | Liouville 函数 $= (-1)^{\Omega(n)}$ | 是（完全积性） |

## 核心理论

### BPSW 素性判定

**BPSW 测试**（Baillie–Pomerance–Selfridge–Wagstaff）将两个独立的概率素数测试串联：一个 base-2 强 Miller–Rabin 测试 + 一个强 Lucas 概率素数测试。自 1980 年提出以来，**没有任何已知合数能通过 BPSW 测试**（$n < 2^{64}$ 范围内已被穷举验证）。

#### Miller–Rabin 测试

**原理**：设 $n - 1 = d \cdot 2^r$（$d$ 为奇数）。对底数 $a$，计算 $a^d \bmod n$。若结果为 $1$ 或在序列 $a^d, a^{2d}, a^{4d}, \ldots, a^{d \cdot 2^{r-1}} \pmod{n}$ 中出现 $-1$，则 $n$ 通过该轮测试。

**判定**：若 $n$ 是素数，则对任意 $1 < a < n$ 都通过。若 $n$ 是合数，则至少 $3/4$ 的底数是证人（witness），即能揭穿 $n$。测试 $k$ 个随机底数后，合数通过全部的概率 $\leq 4^{-k}$。

**确定性**：对 $n < 3.317 \times 10^{24}$，使用固定底数集 $\{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37\}$ 的 Miller–Rabin 测试是**确定性**的。oCAS 使用这 12 个底数作为 `MR_WITNESSES` 常量。

#### 强 Lucas 概率素数测试

**Lucas 序列**：对参数 $(P, Q)$，定义

$$U_0 = 0,\quad U_1 = 1,\quad U_{m+1} = P \cdot U_m - Q \cdot U_{m-1}$$

$$V_0 = 2,\quad V_1 = P,\quad V_{m+1} = P \cdot V_m - Q \cdot V_{m-1}$$

**二进制阶梯（binary ladder）**：高效计算 $(U_k, V_k, Q^k) \bmod n$。利用倍增公式：

$$U_{2m} = U_m \cdot V_m, \quad V_{2m} = V_m^2 - 2Q^m$$

和递增公式：

$$2U_{m+1} = P \cdot U_m + V_m, \quad 2V_{m+1} = D \cdot U_m + P \cdot V_m$$

其中 $D = P^2 - 4Q$。通过模 $n$ 的 half-mod 操作（$n$ 为奇数时 $\text{half\_mod}(x, n) = y$ 满足 $2y \equiv x \pmod{n}$）避免除法。

**Selfridge 参数选择**：选择第一个 $D \in \{5, -7, 9, -11, 13, \ldots\}$ 使得 Jacobi 符号 $\left(\frac{D}{n}\right) = -1$，然后令 $P = 1$，$Q = (1 - D)/4$。

**强 Lucas 测试**：写 $n + 1 = d \cdot 2^r$（$d$ 为奇数），$n$ 通过当且仅当 $U_d \equiv 0 \pmod{n}$ 或存在 $0 \leq i < r$ 使得 $V_{d \cdot 2^i} \equiv 0 \pmod{n}$。

#### BPSW 组合

```
is_prime_bpsw(n):
    1. 若 n < 2，返回 false；若 n 等于 {2,3,5,…,37} 中的素数，返回 true
    2. 若 n 被 2..37 中的某个素数整除，返回 false
    3. 若 n 是完全平方数，返回 false
    4. 强 Miller–Rabin 测试（底数 2）
    5. 若通过，执行强 Lucas PRP 测试（Selfridge 参数）
    6. 两者都通过返回 true
```

oCAS 中 `is_prime()` 使用完整的 12 底数 Miller–Rabin，而 `is_prime_bpsw()` 使用单底数 Miller–Rabin + Lucas 测试。两者互补：`is_prime()` 在 $n < 3.317 \times 10^{24}$ 时确定性；`is_prime_bpsw()` 理论上更强（组合两种不同原理的测试）。

### 整数分解策略

oCAS 的 `factor_integer` 是一个**分层驱动器**（driver），将多种分解算法组合为一个递增策略。

#### 试除法（Trial Division）

首先用 Eratosthenes 筛生成 $\leq 1000$ 的素数列表，逐一试除。这步的复杂度为 $O(\pi(1000)) = O(168)$，代价极小但能快速剥离小因子。

$$n = \prod_{p \leq 1000} p^{e_p} \cdot C$$

余因子 $C$ 的最小素因子 $> 1000$。

#### Pollard rho–Brent 变体

**原理**：基于生日悖论的随机碰撞检测。定义序列 $x_0, x_1, x_2, \ldots$ 满足 $x_{i+1} = x_i^2 + c \pmod{n}$（$c$ 为随机常数）。序列模 $n$ 最终进入循环，而模真因子 $p$ 会更早碰撞。

**Brent 变体**的核心优化：

1. **指数搜索**：用 $2^k$ 步长检测碰撞（Floyd 的龟兔法每步都检测，Brent 只在 2 的幂处检测）
2. **批量 GCD**：积累若干个 $(x_i - y_i)$ 的乘积后一次性取 $\gcd$，摊销 GCD 的高代价
3. **回溯（backtrack）**：当批量 GCD 返回 $n$ 本身时，在该批量内逐个重试

**时间复杂度**：期望 $O(n^{1/4})$ 次模乘。实践中对中等大小因子（< 20 位）非常高效。

#### Pollard $p - 1$ 方法（Stage 1）

**原理**：若 $n$ 有因子 $p$ 且 $p - 1$ 是 $B_1$-光滑的（即 $p - 1$ 的所有素因子幂 $\leq B_1$），则 $p - 1 \mid M$ 其中 $M = \prod_{q \leq B_1} q^{\lfloor \log_q B_1 \rfloor}$。由 Fermat 小定理，对任意 $\gcd(a, p) = 1$：

$$a^M \equiv 1 \pmod{p} \implies p \mid \gcd(a^M - 1, n)$$

**实现细节**：

- 用 Eratosthenes 筛生成 $q \leq B_1$ 的素数
- 逐素数做模幂：$a \leftarrow a^{q^e} \bmod n$（$q^e \leq B_1 < q^{e+1}$）
- 每累积一定步数后检测 $\gcd(a - 1, n)$
- 若 $\gcd = n$（所有因子都 $p-1$ 光滑），方法失败
- 若 $1 < \gcd < n$，找到非平凡因子

**局限**：要求 $p - 1$ 光滑。若 $p - 1$ 有大素因子，方法无效。

#### Williams $p + 1$ 方法（Stage 1）

**原理**：与 $p - 1$ 方法互补。若 $p + 1$ 是 $B_1$-光滑的，利用 Lucas $V$ 序列代替模幂。

**Lucas $V$ 序列**：取 $Q = 1$，随机选择 $P$ 使得 $\left(\frac{P^2 - 4}{n}\right) = -1$（Jacobi 符号）。计算 $V_M \bmod n$，其中 $M = \prod_{q \leq B_1} q^{\lfloor \log_q B_1 \rfloor}$。若 $p + 1$ 是 $B_1$-光滑的，则 $V_M \equiv V_0 = 2 \pmod{p}$，即 $p \mid \gcd(V_M - 2, n)$。

**参数选择**：随机选 $P$，检查 $\text{jacobi}(P^2 - 4, n) = -1$。oCAS 实现中使用 `lucas_uv_mod`（与 BPSW 共享的 Lucas 链代码）。

**局限**：要求 $p + 1$ 光滑。与 $p - 1$ 方法互补（某些 $p$ 满足一个但不满足另一个）。

#### Lenstra 椭圆曲线方法（ECM）

**原理**：将整数分解转化为椭圆曲线群阶的光滑性问题。对随机椭圆曲线 $E/\mathbb{F}_p$，其群阶 $|E(\mathbb{F}_p)|$ 在 Hasse 区间 $[p + 1 - 2\sqrt{p}, p + 1 + 2\sqrt{p}]$ 内均匀分布（Sato–Tate）。ECM 成功的条件是某条曲线的群阶为 $B_1$-光滑的。

**关键优势**：与 rho 和 $p \pm 1$ 方法不同，ECM 不固定群阶区间——每换一条曲线就换一个群阶。这使得 ECM 对**任何大小的因子**都有概率成功，复杂度取决于最小因子 $p$ 的光滑性，而非 $n$ 本身。

**Suyama 参数化**：从随机 $\sigma \notin \{0, 1, 5\}$ 构造 Montgomery 曲线：

$$u = \sigma^2 - 5, \quad v = 4\sigma$$

$$A = \frac{(v - u)^3(3u + v)}{4u^3 v} - 2$$

曲线为 $B y^2 = x^3 + A x^2 + x$，基点为 $(u^3, v^3)$。Suyama 参数化保证基点总是有理的，且在构造过程中若分母的 GCD 暴露出 $n$ 的因子，则立即返回。

**Montgomery 坐标**：使用射影坐标 $(X : Z)$（省略 $Y$ 坐标），曲线方程为 $By^2 = x^3 + Ax^2 + x$，其中 $x = X/Z$。

- **倍点** $[2](X : Z)$：用常量 $a_{24} = (A + 2)/4$ 计算，只需 4 次模乘 + 4 次模加/减
- **差分加法** $P + Q$ 给定 $P - Q$：只需 6 次模乘
- **Montgomery 阶梯**：标量乘 $[k]P$ 通过交替倍点和加法计算，**始终使用基点作为固定差分**

**Stage 1 算法**：计算 $[M]P$ 其中 $M = \prod_{q \leq B_1} q^{\lfloor \log_q B_1 \rfloor}$。在计算过程中取 GCD，若发现非平凡因子则返回。

#### 驱动器 `factor_integer`

`factor_integer` 将上述算法组合为一个递增策略：

```
factor_integer(n):
    1. 试除法：剥离所有 p ≤ 1000 的因子
    2. 对余因子 C（若 C > 1）：
       a. 若 is_prime_bpsw(C)，记录为素因子
       b. 否则调用 find_factor(C) 拆分，递归处理
```

`find_factor(n)` 的升级策略：

1. 先做一次 Pollard rho–Brent（对小因子最廉价，无光滑性要求）
2. 循环（$B_1$ 从 2000 开始）：
   - Pollard $p - 1$（Stage 1，要求某因子 $p - 1$ 为 $B_1$-光滑）
   - Williams $p + 1$（Stage 1，要求某因子 $p + 1$ 为 $B_1$-光滑，与 $p - 1$ 互补）
   - ECM：曲线数 $\approx B_1 / 550$（限制在 10–300 之间）
   - 全部失败后 $B_1 \leftarrow B_1 \times 4$，继续循环

曲线预算约为 $B_1 / 550$ 是 ECM 平滑度试探的标准经验比例，保证每档 $B_1$ 下有足够的随机曲线覆盖 Hasse 区间内的群阶。

### 离散对数

#### Baby-Step Giant-Step（BSGS）

**问题**：给定 $g, h, p$，求 $x$ 使得 $g^x \equiv h \pmod{p}$（$0 \leq x < m$，$m$ 为 $g$ 的阶的上界）。

**算法**：令 $m = \lfloor \sqrt{\text{bound}} \rfloor + 1$。

1. **Baby steps**：计算并存储 $\{(j, g^j \bmod p) : 0 \leq j < m\}$ 到 HashMap
2. **Giant steps**：计算 $g^{-m} \bmod p$，对 $i = 0, 1, \ldots, m-1$，检查 $h \cdot (g^{-m})^i \bmod p$ 是否在表中
3. 若找到匹配 $g^j = h \cdot g^{-im}$，则 $x = im + j$

**复杂度**：时间 $O(\sqrt{m})$ 次模乘，空间 $O(\sqrt{m})$。

oCAS 的 `dlog_bsgs` 使用 `HashMap<Integer, Integer>` 存储 baby steps，用 `(j, g^j)` 键值对实现 $O(1)$ 查找。

#### Pohlig–Hellman 算法

**场景**：$p$ 为素数，$g$ 的阶为 $n = p - 1$，且 $n$ 光滑：$n = q_1^{e_1} \cdots q_k^{e_k}$（$q_i$ 为小素数）。

**核心思想**：将 $x \bmod n$ 的求解分解为 $x \bmod q_i^{e_i}$ 的子问题，然后用 CRT 合并。

**逐素幂分解**：对每个 $q = q_i$，$e = e_i$，求 $x \bmod q^e$：

写 $x = x_0 + x_1 q + x_2 q^2 + \cdots + x_{e-1} q^{e-1}$（$0 \leq x_j < q$，逐"位"恢复）。

1. 令 $g_0 = g^{n/q} \bmod p$，$h_0 = h^{n/q} \bmod p$。则 $g_0^{x_0} = h_0 \pmod{p}$。用 BSGS 求 $x_0 \bmod q$。
2. 更新 $h \leftarrow h \cdot g^{-x_0}$，计算 $h_1 = (h \cdot g^{-x_0})^{n/q^2}$。求 $x_1$。
3. 继续直到 $x_{e-1}$。
4. $x \bmod q^e = x_0 + x_1 q + \cdots + x_{e-1} q^{e-1}$。

**CRT 合并**：对所有 $q_i^{e_i}$ 的结果用 `crt_many` 合并得到 $x \bmod n$。

**复杂度**：每个素幂的 BSGS 花费 $O(\sqrt{q_i})$，总计 $O(\sum e_i \sqrt{q_i})$。当所有 $q_i$ 都小时，这比直接 BSGS 的 $O(\sqrt{n})$ 快得多。

oCAS 的 `dlog_pohlig_hellman` 要求 $p$ 为素数。它使用 `factor_integer` 分解 $p - 1$（复用了整数分解模块），逐素幂执行 BSGS 数字恢复，最后 `crt_many` 合并并验证结果。

### 中国剩余定理

#### 经典 CRT

**定理**：设 $m_1, \ldots, m_k$ 两两互素，给定同余方程组

$$x \equiv r_i \pmod{m_i}, \quad i = 1, \ldots, k$$

则存在唯一解 $x \bmod M$ 其中 $M = m_1 \cdots m_k$。

**构造**：$x = \sum_i r_i \cdot M_i \cdot M_i^{-1}$，其中 $M_i = M / m_i$，$M_i^{-1} = M_i^{-1} \bmod m_i$。

#### 推广：模数不互素

oCAS 的 `crt`（成对合并）和 `crt_many`（多同余合并）**不要求模数两两互素**。

**成对 CRT**：对 $x \equiv r_1 \pmod{m_1}$ 和 $x \equiv r_2 \pmod{m_2}$：

1. 计算 $g = \gcd(m_1, m_2)$
2. 检查 $r_1 \equiv r_2 \pmod{g}$（一致性条件）
3. 若不一致，返回 `None`
4. 否则用扩展 Euclid 算法求解，合并模为 $\text{lcm}(m_1, m_2)$

**多同余合并** `crt_many`：左折叠——从第一个同余开始，逐个与下一个合并：

$$x \equiv r_1 \pmod{m_1} \xrightarrow{\text{合并}} x \equiv R_2 \pmod{\text{lcm}(m_1, m_2)} \xrightarrow{\text{合并}} \cdots$$

最终返回 $(R, M)$ 其中 $M = \text{lcm}(m_1, \ldots, m_k)$，$0 \leq R < M$。若任何一步发现不一致，整体返回 `None`。

**复杂度**：$k$ 次成对合并，每次 $O(\log M)$ 的 GCD 运算。

### 二次剩余

#### Legendre 符号与 Jacobi 符号

**Legendre 符号**：对奇素数 $p$，

$$\left(\frac{a}{p}\right) = \begin{cases} 0 & p \mid a \\ 1 & a \text{ 是模 } p \text{ 二次剩余} \\ -1 & a \text{ 是模 } p \text{ 二次非剩余} \end{cases}$$

**Euler 判据**：$\left(\frac{a}{p}\right) \equiv a^{(p-1)/2} \pmod{p}$。

**Jacobi 符号**：Legendre 符号对合数模的推广。若 $n = p_1^{e_1} \cdots p_k^{e_k}$，则

$$\left(\frac{a}{n}\right) = \prod_i \left(\frac{a}{p_i}\right)^{e_i}$$

**关键性质**：Jacobi 符号可通过**二次互反律**（quadratic reciprocity）快速计算，无需分解 $n$。

#### 二次互反律

对奇素数 $p, q$：

$$\left(\frac{p}{q}\right) \left(\frac{q}{p}\right) = (-1)^{\frac{p-1}{2} \cdot \frac{q-1}{2}}$$

补充规则：

- $\left(\frac{2}{n}\right) = (-1)^{(n^2 - 1)/8}$，即 $n \equiv \pm 1 \pmod 8$ 时为 $1$，$n \equiv \pm 3 \pmod 8$ 时为 $-1$
- $\left(\frac{-1}{n}\right) = (-1)^{(n-1)/2}$，即 $n \equiv 1 \pmod 4$ 时为 $1$

**计算算法**（oCAS 的 `jacobi` 函数）：

```
jacobi(a, n):
    若 n 为偶数或非正，返回 0
    反复：
      1. 2-adic 剥离：提取 a 中所有因子 2，用 mod-8 规则累积符号
      2. 互反律翻转：交换 a ↔ n，根据 a, n 的 mod-4 值调整符号
      3. a = a mod n
      4. 若 a == 0，返回 0（n > 1）或累积符号（n == 1）
```

**复杂度**：$O(\log^2 n)$，与 Euclid 算法同阶。

#### Tonelli–Shanks 算法

求 $x$ 使得 $x^2 \equiv a \pmod{p}$（$p$ 为奇素数，$\left(\frac{a}{p}\right) = 1$）。

**快速路径**：当 $p \equiv 3 \pmod{4}$ 时，

$$x = a^{(p+1)/4} \bmod p$$

因为 $x^2 = a^{(p+1)/2} = a \cdot a^{(p-1)/2} = a \cdot \left(\frac{a}{p}\right) = a$。

**一般情况**（$p \equiv 1 \pmod{4}$）——Tonelli–Shanks 算法：

1. **分解**：写 $p - 1 = Q \cdot 2^S$（$Q$ 为奇数）
2. **找非剩余**：搜索 $z$ 使得 $\left(\frac{z}{p}\right) = -1$
3. **初始化**：
   - $M = S$
   - $c = z^Q \bmod p$（$2^{S-1}$ 阶元素）
   - $t = a^Q \bmod p$
   - $R = a^{(Q+1)/2} \bmod p$
4. **循环**（当 $t \neq 1$）：
   - 找最小 $i > 0$ 使得 $t^{2^i} \equiv 1 \pmod{p}$
   - 更新 $c \leftarrow c^{2^{M-i-1}}$，$t \leftarrow t \cdot c^2$，$R \leftarrow R \cdot c$，$M \leftarrow i$
5. 返回 $R$

**不变量**：$R^2 = a \cdot t \pmod{p}$ 且 $c$ 的阶为 $2^M$。循环终止时 $t = 1$，故 $R^2 \equiv a$。

**复杂度**：$O(S^2 \log p)$，其中 $S = v_2(p - 1)$（$p - 1$ 中 2 的幂次）。平均 $S$ 很小，算法高效。

### 有理重构

**问题**：给定 $a \in \mathbb{Z}$ 和模数 $m$，找到 $n, d \in \mathbb{Z}$ 使得

$$a \cdot d \equiv n \pmod{m}, \quad \gcd(n, d) = 1, \quad 2|n| \cdot |d| < m$$

**应用**：在模 GCD 算法中，通过 CRT 合并多个 $\mathbb{F}_p$ 上的 GCD 系数后，需要将模 $M$ 下的结果"重构"回 $\mathbb{Q}$ 上的有理数。

**Wang/扩展 Euclid 算法**：

使用扩展 Euclid 算法追踪序列 $(r_i, t_i)$：

1. 初始化：$(r_0, r_1) = (m, a)$，$(t_0, t_1) = (0, 1)$
2. 迭代：$q = \lfloor r_0 / r_1 \rfloor$，$(r_0, r_1) \leftarrow (r_1, r_0 - q \cdot r_1)$，$(t_0, t_1) \leftarrow (t_1, t_0 - q \cdot t_1)$
3. **终止条件**：当 $|r_1| \leq \sqrt{m/2}$ 且 $|t_1| \leq \sqrt{m/2}$ 时停止
4. **验证**：令 $n = r_1$，$d = t_1$（取正值），检查 $a \cdot d \equiv n \pmod{m}$ 和 $2|n| \cdot |d| < m$
5. 若验证失败或 $t_1 = 0$，返回 `None`（不存在满足条件的有理重构）

**唯一性定理**：当 $2|n| \cdot |d| < m$ 时，满足条件的 $(n, d)$ 是唯一的。这个条件保证了重构结果是"最简"的有理表示。

**复杂度**：$O(\log m)$ 的 GCD 运算步数，与 Euclid 算法同阶。

## 在 oCAS 中的实现

oCAS 的数论算法实现在 `ocas-domain` crate 的 `number_theory` 模块和 `ocas-poly` 的有理重构模块中。

### 模块结构

```
ocas-domain/src/
├── number_theory.rs          ← 素性测试、模逆、CRT、Legendre/Jacobi、mod_sqrt
└── number_theory/
    ├── primes.rs             ← BPSW、Lucas 序列（lucas_uv_mod、strong_lucas_prp）
    ├── factor.rs             ← 整数分解驱动器（factor_integer）及各方法
    ├── dlog.rs               ← BSGS 和 Pohlig–Hellman
    ├── crt.rs                ← 多同余 CRT（crt_many）
    └── functions.rs          ← φ、μ、τ、σ、λ

ocas-poly/src/
└── rational_reconstruction.rs ← 有理重构（扩展 Euclid 方法）
```

### 素性判定

| 函数 | 位置 | 说明 |
|---|---|---|
| `is_prime(n)` | `number_theory.rs` | 12 底数 Miller–Rabin，$n < 3.317 \times 10^{24}$ 确定性 |
| `is_prime_bpsw(n)` | `primes.rs` | base-2 MR + 强 Lucas PRP |
| `is_prime_u64(n)` | `primes.rs` | `u64` 专用，委托给 `is_prime` |
| `next_prime(n)` | `number_theory.rs` | 从 $n + 1$ 起逐个奇数测试 |
| `primes_from(n)` | `number_theory.rs` | 素数迭代器 |

Lucas 序列计算由 `lucas_uv_mod` 实现，使用二进制阶梯法，返回 $(U_k, V_k, Q^k) \bmod n$。该函数同时被 BPSW 测试和 Williams $p + 1$ 分解方法共用。

### 整数分解

| 函数 | 位置 | 说明 |
|---|---|---|
| `factor_integer(n)` | `factor.rs` | 入口：试除 + 递归拆分 |
| `factor_integer_with_rng(n, rng)` | `factor.rs` | 带显式 RNG 的版本 |
| `factor_trial(n, limit)` | `factor.rs` | 试除法，返回 `(因子, 余因子)` |
| `pollard_rho_brent(n, rng)` | `factor.rs` | Brent 变体，带重试 |
| `pollard_pm1(n, b1, rng)` | `factor.rs` | Pollard $p - 1$，Stage 1 |
| `williams_pp1(n, b1, rng)` | `factor.rs` | Williams $p + 1$，Stage 1 |
| `ecm(n, b1, max_curves, rng)` | `factor.rs` | Lenstra ECM，Suyama 参数化 |

内部实现细节：

- `primes_up_to(limit)`: Eratosthenes 筛，为试除法和光滑界计算生成素数表
- `prime_power_le(q, bound)`: 计算 $q^e \leq \text{bound}$，用于构造光滑指数 $M$
- `ProjPoint { x, z }`: Montgomery 射影坐标点
- `ecm_double`、`ecm_add`、`ecm_mul`: Montgomery 曲线上的群运算
- `suyama_curve(sigma, n)`: Suyama 参数化，返回 `Suyama` 枚举（`Curve`/`Factor`/`Degenerate`）
- `find_factor(n, rng)`: 升级策略驱动器（rho → p−1 → p+1 → ECM，光滑界递增）

### 离散对数

| 函数 | 位置 | 说明 |
|---|---|---|
| `dlog_bsgs(base, target, modulus)` | `dlog.rs` | Baby-step giant-step |
| `dlog_pohlig_hellman(base, target, p)` | `dlog.rs` | Pohlig–Hellman（需 $p$ 为素数）|

内部：

- `bsgs_bounded(base, target, modulus, order_bound)`: BSGS 核心，给定阶上界搜索
- `dlog_bsgs` 设置上界为 `modulus - 1` 并委托给 `bsgs_bounded`
- `dlog_pohlig_hellman` 使用 `factor_integer` 分解 $p - 1$，逐素幂 BSGS 恢复数字，`crt_many` 合并

### 中国剩余定理

| 函数 | 位置 | 说明 |
|---|---|---|
| `crt(r1, m1, r2, m2)` | `number_theory.rs` | 成对合并，模数无需互素 |
| `crt_many(congruences)` | `crt.rs` | 多同余左折叠合并 |

`crt_many` 从第一个同余开始，逐个调用 `crt` 合并。任何一步不一致则整体返回 `None`。

### 二次剩余

| 函数 | 位置 | 说明 |
|---|---|---|
| `legendre(a, p)` | `number_theory.rs` | Legendre 符号，委托给 `jacobi` |
| `jacobi(a, n)` | `number_theory.rs` | Jacobi 符号（二次互反律 + 2-adic 剥离） |
| `mod_sqrt(a, p)` | `number_theory.rs` | Tonelli–Shanks，含 $p \equiv 3 \pmod{4}$ 快速路径 |

`mod_sqrt` 实现：

1. 先检查 $\left(\frac{a}{p}\right) = 1$（否则返回 `None`）
2. 若 $p \equiv 3 \pmod{4}$，走快速路径 $x = a^{(p+1)/4} \bmod p$
3. 否则执行完整 Tonelli–Shanks 算法：找非剩余 $z$，分解 $p - 1 = Q \cdot 2^S$，维护 $(c, t, r, m)$ 不变量循环

### 辅助数论函数

| 函数 | 位置 | 说明 |
|---|---|---|
| `euler_phi(n)` | `functions.rs` | $\varphi(n) = |n| \prod(1 - 1/p)$，基于 `factor_integer` |
| `moebius_mu(n)` | `functions.rs` | $\mu(n)$，基于 `factor_integer` |
| `divisor_tau(n)` | `functions.rs` | $\tau(n) = \prod(e_i + 1)$ |
| `divisor_sigma(n, k)` | `functions.rs` | $\sigma_k(n) = \prod \frac{p^{k(e+1)} - 1}{p^k - 1}$ |
| `liouville_lambda(n)` | `functions.rs` | $\lambda(n) = (-1)^{\Omega(n)}$ |

所有积性函数都先通过 `factor_integer` 获取素因子分解，再用积性公式计算。

### 有理重构

| 函数 | 位置 | 说明 |
|---|---|---|
| `rational_reconstruction(a, m)` | `ocas-poly/src/rational_reconstruction.rs` | Wang/扩展 Euclid 方法 |

内部使用 `integer_sqrt` 计算 $\lfloor\sqrt{m/2}\rfloor$ 作为终止界的阈值。

### 基础工具函数

| 函数 | 位置 | 说明 |
|---|---|---|
| `mod_inv(a, m)` | `number_theory.rs` | 扩展 Euclid 求模逆 |
| `extended_gcd(a, b)` | `number_theory.rs` | 返回 $(g, x, y)$ 满足 $ax + by = g$ |
| `symmetric_mod(a, m)` | `number_theory.rs` | 约化到 $(-m/2, m/2]$ |

## 参考文献

1. **Shoup, V.** *A Computational Introduction to Number Theory and Algebra.* Cambridge University Press, 2nd edition, 2009.
   - Ch.10: Primality testing (Miller–Rabin, randomized algorithms)
   - Ch.11: Finding discrete logarithms (BSGS, Pohlig–Hellman)
   - Ch.19: Factoring integers (Pollard rho, ECM)

2. **Crandall, R. & Pomerance, C.** *Prime Numbers: A Computational Perspective.* Springer, 2nd edition, 2005.
   - Ch.3: Recognizing primes and composites (BPSW, Lucas tests)
   - Ch.5: Exponential factoring algorithms (ECM, Pollard rho)
   - Ch.6: Subexponential factoring (Pollard $p-1$, Williams $p+1$)
   - Ch.7: Modern discrete logarithm algorithms

3. **Brent, R. P.** "An improved Monte Carlo factorization algorithm." *BIT Numerical Mathematics*, 20(2):176–184, 1980.
   - Brent 变体的 Pollard rho 算法

4. **Williams, H. C.** "A p+1 method of factoring." *Mathematics of Computation*, 39(159):225–234, 1982.
   - Williams $p + 1$ 方法

5. **Lenstra, H. W. Jr.** "Factoring integers with elliptic curves." *Annals of Mathematics*, 126(3):649–673, 1987.
   - ECM 方法

6. **Montgomery, P. L.** "Speeding the Pollard and elliptic curve methods of factorization." *Mathematics of Computation*, 48(177):243–264, 1987.
   - Montgomery 曲线坐标与阶梯优化

7. **Shanks, D. & Tonelli, R.** See Cohen, H. *A Course in Computational Algebraic Number Theory*, Algorithm 1.5.1.
   - Tonelli–Shanks 模平方根算法

8. **Baillie, R. & Wagstaff, S. S.** "Lucas pseudoprimes." *Mathematics of Computation*, 35(152):1391–1417, 1980.
   - BPSW 测试的 Lucas 部分

9. **Pohlig, S. & Hellman, M.** "An improved algorithm for computing logarithms over GF(p) and its cryptographic significance." *IEEE Transactions on Information Theory*, 24(1):106–110, 1978.
   - Pohlig–Hellman 离散对数算法

10. **Wang, P. S.** "A p-adic algorithm for univariate partial fractions." *Proceedings of SYMSAC '81*, 212–217, 1981.
    - 有理重构的扩展 Euclid 方法
