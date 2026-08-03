# 高阶：蒙特卡洛积分

## 前提知识

- [线性代数](./linear-algebra.md) — 矩阵运算（用于理解高维采样）
- 基础概率论 — 随机变量、期望、方差
- 基础微积分 — 定积分的定义与变量替换

建议具备以下概率论基础：

| 概念 | 说明 |
|------|------|
| 期望 $\mathbb{E}[X]$ | 随机变量的加权平均值 |
| 方差 $\text{Var}(X) = \mathbb{E}[(X - \mu)^2]$ | 随机变量偏离均值的度量 |
| 独立同分布 (i.i.d.) | 多个随机变量相互独立且服从同一分布 |
| 大数定律 | 样本均值随样本量增大收敛于期望 |
| 中心极限定理 | 样本均值近似服从正态分布 |

---

## 基础概念

### 从定积分到期望

考虑定积分

$$I = \int_\Omega f(\mathbf{x})\, d\mathbf{x}$$

其中 $\Omega \subseteq \mathbb{R}^d$ 是积分域。蒙特卡洛方法的核心思想是将积分**重写为期望**。

令 $\mathbf{X}$ 为在 $\Omega$ 上服从均匀分布的随机变量，其概率密度函数为

$$p(\mathbf{x}) = \frac{1}{V}, \qquad V = \text{vol}(\Omega)$$

则积分变为

$$I = V \cdot \mathbb{E}[f(\mathbf{X})]$$

### 基本估计量

**蒙特卡洛基本估计量**：从 $p(\mathbf{x})$ 中抽取 $N$ 个独立样本 $\mathbf{x}_1, \dots, \mathbf{x}_N$，定义

$$\hat{I} = \frac{V}{N} \sum_{i=1}^{N} f(\mathbf{x}_i)$$

**性质**：

1. **无偏性**：$\mathbb{E}[\hat{I}] = I$（精确等于真值的期望）。
2. **方差**：$\text{Var}(\hat{I}) = \frac{V^2}{N} \text{Var}(f(\mathbf{X}))$。
3. **收敛速率**：$\hat{I}$ 的标准误差为 $\sigma_{\hat{I}} \propto 1/\sqrt{N}$，**与维度 $d$ 无关**。

最后一条性质是蒙特卡洛方法的核心优势：传统数值积分（梯形法则、Simpson 法则）的误差随维度指数增长（"维度灾难"），而蒙特卡洛方法的收敛速率始终为 $O(N^{-1/2})$。

### 大数定律与中心极限定理

**大数定律**（Kolmogorov 强大数定律）. 若 $X_1, X_2, \dots$ 是独立同分布的随机变量，$\mathbb{E}[|X_1|] < \infty$，则

$$\bar{X}_n = \frac{1}{n}\sum_{i=1}^{n} X_i \;\xrightarrow{\text{a.s.}}\; \mathbb{E}[X_1]$$

应用到蒙特卡洛：$\hat{I}/V \to \mathbb{E}[f(\mathbf{X})]$ 几乎必然成立，因此 $\hat{I} \to I$。

**中心极限定理**（CLT）. 在上述条件下，若 $\text{Var}(X_1) = \sigma^2 < \infty$，则

$$\frac{\bar{X}_n - \mu}{\sigma / \sqrt{n}} \;\xrightarrow{d}\; \mathcal{N}(0, 1)$$

应用到蒙特卡洛：对足够大的 $N$，$\hat{I}$ 近似服从正态分布：

$$\hat{I} \;\sim\; \mathcal{N}\!\left(I,\; \frac{V^2 \sigma_f^2}{N}\right)$$

其中 $\sigma_f^2 = \text{Var}(f(\mathbf{X}))$。这给出了**置信区间**的基础。

### 方差估计与置信区间

实际中 $\sigma_f^2$ 未知，用样本方差估计：

$$\hat{\sigma}_f^2 = \frac{1}{N-1} \sum_{i=1}^{N} \bigl(f(\mathbf{x}_i) - \bar{f}\bigr)^2, \qquad \bar{f} = \frac{1}{N}\sum_{i=1}^{N} f(\mathbf{x}_i)$$

$\hat{I}$ 的标准误差估计为

$$\hat{\sigma}_{\hat{I}} = \frac{V}{\sqrt{N}} \hat{\sigma}_f$$

$95\%$ 置信区间为 $\hat{I} \pm 1.96\, \hat{\sigma}_{\hat{I}}$。

**加权采样的推广**：当采样密度 $p(\mathbf{x})$ 不是均匀分布时（重要性采样），每个样本 $i$ 带有权重 $w_i = 1/p(\mathbf{x}_i)$，则

$$\hat{I} = \frac{1}{N} \sum_{i=1}^{N} w_i \cdot f(\mathbf{x}_i)$$

方差估计变为

$$\hat{\sigma}^2 = \frac{1}{N} \sum_{i=1}^{N} w_i^2 f(\mathbf{x}_i)^2 - \left(\frac{1}{N}\sum_{i=1}^{N} w_i f(\mathbf{x}_i)\right)^2$$

oCAS 的 `StatisticsAccumulator` 使用与之类似的加权方差形式（见下文『统计累加器』一节）。

---

## 核心理论

### 重要性采样

蒙特卡洛估计量的方差取决于被积函数在采样分布下的"波动"程度。**重要性采样**（importance sampling）通过选择一个非均匀的采样分布 $p(\mathbf{x})$ 来降低方差。

**原理**：将积分重写为

$$I = \int_\Omega \frac{f(\mathbf{x})}{p(\mathbf{x})} \cdot p(\mathbf{x})\, d\mathbf{x} = \mathbb{E}_p\!\left[\frac{f(\mathbf{X})}{p(\mathbf{X})}\right]$$

估计量为

$$\hat{I}_p = \frac{1}{N} \sum_{i=1}^{N} \frac{f(\mathbf{x}_i)}{p(\mathbf{x}_i)}, \qquad \mathbf{x}_i \sim p$$

**最优采样分布**：方差最小化的 $p^*$ 为

$$p^*(\mathbf{x}) = \frac{|f(\mathbf{x})|}{\int_\Omega |f(\mathbf{y})|\, d\mathbf{y}}$$

当 $f$ 不变号时，$w \cdot f = \int |f| \cdot \mathrm{sgn}(f)$ 是常数，此时方差为零；$f$ 变号时方差也达到最小。但 $p^*$ 本身需要知道积分值——这是一个"先有鸡还是先有蛋"的问题。实际中，我们用**自适应网格**来近似 $p^*$。

**关键约束**：$p(\mathbf{x})$ 必须在 $f(\mathbf{x}) \neq 0$ 的所有点上严格为正。若 $p$ 在 $f$ 非零处为零，估计量无界，方差无穷大。

### Vegas 自适应网格

Vegas 算法（Lepage, 1978）是最广泛使用的自适应蒙特卡洛积分方法。其核心思想是用一个**分段常数**的采样密度来近似最优分布 $p^*$。

#### 乘积网格结构

Vegas 在 $d$ 维单位超立方体 $[0,1]^d$ 上工作（物理域通过线性变量替换映射到超立方体）。采样密度采用**乘积形式**：

$$p(\mathbf{x}) = \prod_{k=1}^{d} p_k(x_k)$$

每个 $p_k$ 是第 $k$ 维上的独立一维分段常数密度。这意味着每维的网格独立维护——$d$ 维问题被分解为 $d$ 个一维问题。

**第 $k$ 维的网格**：将 $[0,1]$ 分为 $M$ 个 bin，边界为

$$0 = b_0 < b_1 < \cdots < b_M = 1$$

第 $j$ 个 bin 的宽度为 $\Delta_j = b_{j+1} - b_j$。在该 bin 上 $p_k(x) = 1/(M \cdot \Delta_j)$（常数），因为每个 bin 被等概率选中。

#### 采样过程

对每个样本点 $\mathbf{x} = (x_1, \dots, x_d)$：

1. **选 bin**：对第 $k$ 维，从 $\{0, 1, \dots, M-1\}$ 中均匀随机选取 bin 编号 $j_k$。
2. **选点**：在 $[b_{j_k}, b_{j_k+1}]$ 内均匀随机选取 $x_k$。
3. **计算 Jacobian（权重）**：该点的概率密度为 $p_k(x_k) = 1/(M \cdot \Delta_{j_k})$，因此

$$w_k = \frac{1}{p_k(x_k)} = M \cdot \Delta_{j_k}$$

总权重为

$$w = \prod_{k=1}^{d} w_k = \prod_{k=1}^{d} (M \cdot \Delta_{j_k})$$

4. **累加**：将 $w \cdot f(\mathbf{x})$ 加入统计累加器。

**权重的物理意义**：$w = 1/p(\mathbf{x})$ 是采样密度的倒数——稀疏区域（$\Delta_j$ 大）的样本获得更大权重，密集区域（$\Delta_j$ 小）的样本获得更小权重，从而补偿了非均匀采样。

#### 网格自适应：训练信号

每次迭代中，Vegas 记录每个 bin 的"重要性"。对第 $k$ 维的第 $j$ 个 bin，训练信号为

$$d_j = \sum_{\text{samples in bin } j} w \cdot f^2(\mathbf{x})$$

**直觉**：$f^2$ 度量该 bin 区域内被积函数的"能量"（方差贡献），$w$ 补偿非均匀采样。高 $d_j$ 的 bin 对积分误差贡献大，需要更细的网格。

#### 网格更新：平滑累积弧长重分布

有了训练信号 $\{d_0, d_1, \dots, d_{M-1}\}$，Vegas 通过以下步骤更新网格边界：

**步骤 1：3-bin 平滑**。对每个 bin $j$，计算平滑值：

$$\tilde{d}_j = \frac{d_{j-1} + d_j + d_{j+1}}{3}$$

边界 bin（$j = 0$ 或 $j = M-1$）缺失的邻居用 0 代替。平滑的目的是抑制单 bin 的噪声，使网格变化更稳定。

**步骤 2：均值归一化**。计算平均值 $\bar{d} = \frac{1}{M}\sum_j \tilde{d}_j$，令

$$\hat{d}_j = \frac{\tilde{d}_j}{\bar{d}}$$

归一化后 $\hat{d}_j = 1$ 表示"平均重要性"，$\hat{d}_j > 1$ 表示该 bin 需要更细的网格。

**步骤 3：阻尼**。引入学习率 $\alpha > 0$（oCAS 默认 $\alpha = 1.5$），对归一化值进行幂变换：

$$\hat{d}_j \;\leftarrow\; \hat{d}_j^{1/\alpha}$$

当 $\alpha = 1$ 时无阻尼（$\hat{d}_j^{1/1} = \hat{d}_j$）；$\alpha > 1$ 时幂次 $< 1$，压缩极端值，使网格更新更保守。这防止了网格对单次迭代的噪声过度反应。

**步骤 4：累积弧长与等距重分布**。计算累积和：

$$S_j = \sum_{i=0}^{j-1} \hat{d}_i, \qquad j = 0, 1, \dots, M$$

其中 $S_0 = 0$，$S_M = \sum_{i=0}^{M-1} \hat{d}_i$。将 $S$ 视为一条"弧长"曲线，新的 bin 边界在弧长上等距分布：

$$b'_j = \text{在 } [0,1] \text{ 上对应弧长 } S = j \cdot S_M / M \text{ 的位置}$$

具体地，对 $j = 1, \dots, M-1$，找到 $S$ 中满足 $S_i \leq j \cdot S_M / M < S_{i+1}$ 的位置 $i$，然后线性插值：

$$b'_j = \frac{i + \frac{j \cdot S_M / M - S_i}{S_{i+1} - S_i}}{M}$$

固定 $b'_0 = 0$，$b'_M = 1$。

**步骤 5：单调性钳制**。由于浮点精度和数值噪声，$b'_j$ 可能不严格递增。对结果进行单调性修复：

$$b'_j \;\leftarrow\; \max(b'_j,\; b'_{j-1})$$

**直觉**：弧长重分布的几何意义是——如果我们在 $(x, S(x))$ 平面上画出累积重要性曲线，"等弧长分点"恰好使得每个新 bin 在曲线上的弧长相等，即每个新 bin 承载近似相等的积分贡献。

### 分层采样

**分层采样**（stratified sampling）是另一种方差缩减技术。其思想是将积分域划分为不相交的子域（"层"），在每层内独立采样：

$$I = \sum_{j=1}^{M} I_j, \qquad I_j = \int_{\Omega_j} f(\mathbf{x})\, d\mathbf{x}$$

每层的估计量为 $\hat{I}_j = V_j \cdot \bar{f}_j$（$V_j$ 是层体积，$\bar{f}_j$ 是层内样本均值），总估计量为 $\hat{I} = \sum_j \hat{I}_j$。

**方差分析**：

$$\text{Var}(\hat{I}) = \sum_{j=1}^{M} \frac{V_j^2}{N_j} \sigma_j^2$$

其中 $N_j$ 是第 $j$ 层的样本数，$\sigma_j^2$ 是第 $j$ 层内 $f$ 的方差。若各层方差不等，最优分配（Neyman 分配）为

$$N_j \propto V_j \sigma_j$$

即高方差层分配更多样本。

**Vegas 与分层采样的关系**：Vegas 的乘积网格本质上是一种**自适应分层采样**——bin 边界将积分域分层，每层内均匀采样。区别在于 Vegas 的层边界是迭代自适应的，而非事先固定的。

### 多轮迭代的逆方差加权组合

Vegas 运行多轮迭代（每轮 $N$ 个样本），每轮给出一个估计 $\hat{I}_r$ 及其标准误差 $\hat{\sigma}_r$。如何组合多轮结果？

**逆方差加权**（inverse-variance weighting）：令 $w_r = 1/\hat{\sigma}_r^2$，组合估计为

$$\hat{I}_{\text{combined}} = \frac{\sum_r w_r \hat{I}_r}{\sum_r w_r}$$

组合误差为

$$\hat{\sigma}_{\text{combined}} = \frac{1}{\sqrt{\sum_r w_r}}$$

**最优性**：当各轮估计相互独立时，逆方差加权是**方差最小**的无偏组合方式（Gauss–Markov 定理的直接推论）。

**在线逆方差加权累加**：为避免存储所有轮次的估计值，`StatisticsAccumulator` 使用在线更新（加权均值 + 增量 $\chi^2$，并非 Welford 算法本身）。维护状态变量：

- `integral`：当前组合估计
- `error`：当前组合标准误差
- `chi_square`：跨轮次 $\chi^2$ 统计量

当第 $r$ 轮完成（估计 $\mu_r$，误差 $\sigma_r$）时，更新规则为：

$$w_{\text{prev}} = 1/\sigma_{\text{prev}}^2, \quad w_r = 1/\sigma_r^2, \quad w_{\text{new}} = w_{\text{prev}} + w_r$$

$$\hat{I}_{\text{new}} = \frac{w_{\text{prev}} \cdot \hat{I}_{\text{prev}} + w_r \cdot \mu_r}{w_{\text{new}}}$$

$$\sigma_{\text{new}} = 1/\sqrt{w_{\text{new}}}$$

**实现细节**：为防止零方差轮次导致权重溢出，oCAS 将误差钳制在 $\sigma_r \geq 10^{-150}$（$10^{-300}$ 仍可表示为 `f64`）。

### $\chi^2$ 诊断

**目的**：检验各轮次的估计是否一致——若各轮的 $\hat{I}_r$ 差异远超其误差条，说明方差估计偏低或积分器未收敛。

**定义**：

$$\chi^2 = \sum_{r=1}^{R} w_r (\hat{I}_r - \hat{I}_{\text{combined}})^2$$

其中 $w_r = 1/\hat{\sigma}_r^2$。

**解释**：若各轮独立且方差估计准确，$\chi^2$ 服从自由度为 $R - 1$ 的卡方分布。经验准则：

| $\chi^2 / (R-1)$ | 含义 |
|:---:|------|
| $\approx 1$ | 正常——各轮一致 |
| $\gg 1$ | 方差被低估，或积分器未收敛 |
| $\ll 1$ | 方差被高估（罕见） |

oCAS 的 `StatisticsAccumulator` 在每轮结束时在线更新 $\chi^2$：

$$\chi^2_{\text{new}} = \chi^2_{\text{prev}} + w_{\text{prev}} (\hat{I}_{\text{prev}} - \hat{I}_{\text{new}})^2 + w_r (\mu_r - \hat{I}_{\text{new}})^2$$

注意这是增量式更新——无需存储历史轮次的估计值。

---

## 在 oCAS 中的实现

oCAS 的 Vegas 实现位于 `ocas-eval/src/numeric/vegas.rs`，配套的统计累加器在 `ocas-eval/src/numeric/statistics.rs`。

### 数据结构

```
Vegas
├── opts: VegasOptions        // 调参选项
├── axes: Vec<GridAxis>       // 每维独立的网格
│   ├── boundaries: Vec<f64>  // M+1 个 bin 边界 ∈ [0,1]
│   └── bin_accum: Vec<f64>   // M 个 bin 的训练信号累加器
└── accumulator: StatisticsAccumulator  // 跨轮次逆方差加权
    ├── sum_w, sum_wf, sum_wf2  // 当前轮次的加权累加
    ├── integral, error          // 组合估计与误差
    ├── chi_square               // χ² 诊断
    └── iterations               // 已完成轮次数
```

### 默认参数

| 参数 | 默认值 | 含义 |
|------|:------:|------|
| `n_bins` | 64 | 每维 bin 数 $M$ |
| `n_samples` | 10,000 | 每轮采样数 $N$ |
| `iterations` | 10 | 自适应轮次数 $R$ |
| `learning_rate` | 1.5 | 网格更新阻尼参数 $\alpha$ |
| `seed` | `0x0C45` | 确定性 RNG 种子 |

### 采样过程的实现

每轮迭代的伪代码：

```rust
for _ in 0..n_samples {
    // 1. 对每维独立采样
    let mut x = Vec::with_capacity(n_dims);
    let mut jac = 1.0;
    for axis in axes.iter_mut() {
        // 均匀选 bin
        let b = rng.random_range(0..n_bins);
        let lo = boundaries[b];
        let hi = boundaries[b + 1];
        // bin 内均匀选点
        let u = rng.random::<f64>();
        let xi = lo + (hi - lo) * u;
        // Jacobian = M × bin_width = 1/pdf
        let wi = (hi - lo) * n_bins as f64;
        x.push(xi);
        jac *= wi;
    }
    // 2. 求值并累加
    let fx = f(&x);
    accumulator.add_sample(jac, fx);
    // 3. 记录训练信号
    for (i, xi) in x.iter().enumerate() {
        axes[i].add_training(*xi, jac, fx * fx);
    }
}
```

**关键实现细节**：

- **RNG 选择**：使用 `Xoshiro256PlusPlus`（而非 `thread_rng`），通过固定种子保证结果可重现。
- **乘积 Jacobian**：总权重 `jac = ∏ w_k`——各维独立采样，总密度为各维密度之积。
- **训练信号**：`add_training` 使用二分查找定位 bin（$O(\log M)$），将 `weight × f²` 累加到对应 bin。

### 网格更新的实现

`GridAxis::update(learning_rate)` 的实现步骤：

```rust
fn update(&mut self, learning_rate: f64) {
    let n = self.bin_accum.len();
    let total: f64 = self.bin_accum.iter().sum();
    if total <= 0.0 { return; }

    // 1. 3-bin 平滑
    let avg = total / n as f64;
    let mut d = vec![0.0; n];
    for i in 0..n {
        let prev = if i > 0 { bin_accum[i-1] } else { 0.0 };
        let next = if i+1 < n { bin_accum[i+1] } else { 0.0 };
        d[i] = (prev + bin_accum[i] + next) / 3.0 / avg;
    }

    // 2. 阻尼：d[i] = d[i]^(1/learning_rate)
    if learning_rate != 1.0 {
        for v in d.iter_mut() {
            *v = v.max(1e-30).powf(1.0 / learning_rate);
        }
    }

    // 3. 累积弧长
    let mut cum = vec![0.0; n+1];
    for i in 0..n { cum[i+1] = cum[i] + d[i]; }

    // 4. 等距重分布
    let mut new_boundaries = vec![0.0; n+1];
    new_boundaries[0] = 0.0;
    new_boundaries[n] = 1.0;
    let mut j = 0;
    for i in 1..n {
        let target = i as f64 / n as f64 * cum[n];
        while j < n && cum[j+1] < target { j += 1; }
        let frac = (target - cum[j]) / (cum[j+1] - cum[j]);
        new_boundaries[i] = (j as f64 + frac) / n as f64;
    }

    // 5. 单调性钳制
    for i in 1..=n {
        if new_boundaries[i] < new_boundaries[i-1] {
            new_boundaries[i] = new_boundaries[i-1];
        }
    }
    new_boundaries[n] = 1.0;
    self.boundaries = new_boundaries;
    self.bin_accum.fill(0.0);  // 重置训练信号
}
```

**实现中的数值保护**：

- 归一化后的值在幂变换前先取 `max(1e-30)`，保证 $d_i > 0$，使弧长累积严格递增（防止后续插值除零/退化）。
- 当 `learning_rate = 1.0` 时幂变换（连同 1e-30 钳制）被跳过，此时 `cum[j+1] - cum[j]` 仍可能为零（某 bin 完全平坦），`frac` 取 0。
- 最终钳制保证边界严格非递减，首尾固定为 0 和 1。

### 统计累加器

`StatisticsAccumulator` 在每轮中用三个变量在线累加：

$$S_w = \sum_i w_i, \qquad S_{wf} = \sum_i w_i f_i, \qquad S_{wf^2} = \sum_i w_i f_i^2$$

轮次结束时，估计量为

$$\mu = \frac{S_{wf}}{S_w}, \qquad \sigma^2 = \frac{S_{wf^2}}{S_w} - \mu^2$$

注意这不是标准的样本方差公式——这是**加权方差**公式，适用于每个样本有不同权重的情况（Vegas 的 Jacobian 权重）。

### 一维便捷函数

`integrate_1d(f, a, b, opts)` 内部做线性变量替换 $x = a + u \cdot (b-a)$：

```rust
let width = b - a;
let wrapped = |u: &[f64]| f(a + u[0] * width) * width;
let mut vegas = Vegas::new(1, opts);
vegas.integrate(&wrapped)
```

Jacobian 因子 $(b-a)$ 直接乘入被积函数，用户传入的 $f(x)$ 接收物理坐标而非超立方体坐标。

---

## 进阶话题

### 何时使用 Vegas vs 符号积分

| 场景 | 推荐方法 |
|------|----------|
| 被积函数有初等原函数 | 符号积分（`integrate`）——精确结果 |
| 被积函数含特殊函数（Bessel、误差函数等） | 数值积分（Vegas） |
| 高维积分（$d > 3$） | 蒙特卡洛（Vegas 或其他） |
| 被积函数有奇点或边界层 | Vegas + 增加 `n_bins` 和 `iterations` |
| 需要可微分的积分值 | 自动微分（见 [自动微分](./autodiff.md)） |

### Vegas 的局限性

1. **乘积网格假设**：Vegas 假设最优采样密度可分离为 $\prod_k p_k(x_k)$。当被积函数有强变量间相关性（如 $f(x,y) = \delta(x - y)$）时，乘积网格无法有效近似，收敛变慢。

2. **分段常数近似**：网格密度是分段常数的——在每个 bin 内采样均匀。若被积函数在某个 bin 内有剧烈变化（如窄峰），该 bin 的方差估计可能不准确。

3. **自适应偏差**：早期迭代使用尚未收敛的网格，其样本不完全独立。Vegas 通过逆方差加权缓解此问题——收敛后的迭代权重更大。

4. **确定性种子**：oCAS 使用固定 RNG 种子（默认 `0x0C45`）保证可重现性。这在测试和基准中有用，但生产环境中应考虑使用随机种子以避免系统性偏差。

### 改进方向

- **Vegas+**（Lepage, 2020）：在网格更新中引入"反射"策略，处理对称被积函数。
- **Suave**（Hahn, 2005）：结合 Vegas 自适应网格与分层采样的优势。
- **Divonne**（Friedman & Harris, 1996）：使用分区树和启发式细分策略。
- **Cuhre**：纯确定性求积（非蒙特卡洛），适用于低维光滑被积函数。

这些方法可通过实现 `Integrator` trait 扩展到 oCAS 中。

---

## 参考文献

1. **Lepage, G. P.** "A new algorithm for adaptive multidimensional integration." *Journal of Computational Physics*, 27(2):192–203, 1978. — Vegas 算法的原始论文。

2. **James, F.** "Monte Carlo theory and practice." *Reports on Progress in Physics*, 43(9):1145–1189, 1980. — 蒙特卡洛方法的综述，包含 Vegas 的详细描述。

3. **Lepage, G. P.** "Vegas: An adaptive multi-dimensional integration program." *Cornell preprint CLNS 80-447*, 1980. — Vegas 的实现细节和使用指南。

4. **Hahn, T.** "CUBA—a library for multidimensional numerical integration." *Computer Physics Communications*, 168(2):78–95, 2005. — Vegas、Suave、Divonne、Cuhre 的比较。

5. **Neal, R. M.** "Annealed importance sampling." *Statistics and Computing*, 11(2):125–139, 2001. — 重要性采样的理论基础。

6. **Kahn, H. & Marshall, A. W.** "Methods of reducing sample size in Monte Carlo computations." *Journal of the Operations Research Society of America*, 1(5):263–278, 1953. — 分层采样和重要性采样的早期工作。
