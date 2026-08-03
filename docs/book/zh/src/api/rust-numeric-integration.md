# Rust API 参考：数值积分

本章涵盖 `ocas-eval` crate 中的数值积分功能，核心是 **Vegas 自适应蒙特卡洛积分器**（Lepage 1978）。所有接口通过 `ocas_eval::numeric` 模块导出。

**模块概览**

```rust
use ocas_eval::numeric::{
    Vegas,                // 自适应蒙特卡洛积分器
    VegasOptions,         // 配置参数
    IntegrateResult,      // 积分结果
    Integrator,           // 积分器 trait
    integrate_1d,         // 一维便捷函数
    StatisticsAccumulator, // 逆方差加权统计累加器
};
```

---

## IntegrateResult

**签名**：

```rust
#[derive(Debug, Clone, Copy)]
pub struct IntegrateResult {
    pub integral: f64,
    pub error: f64,
}
```

**功能**：数值积分的返回结果，包含积分估计值和标准误差。

**字段**：

| 字段 | 类型 | 说明 |
|---|---|---|
| `integral` | `f64` | 积分的最佳估计值 |
| `error` | `f64` | `integral` 的标准误差估计 |

**示例**：

```rust
use ocas_eval::numeric::integrate_1d;

let r = integrate_1d(|x| x * x, 0.0, 1.0, Default::default());
println!("∫₀¹ x² dx ≈ {:.6} ± {:.6}", r.integral, r.error);
// 输出：∫₀¹ x² dx ≈ 0.3334 ± 0.0003 （精确值 1/3）
```

---

## Integrator

**签名**：

```rust
pub trait Integrator {
    fn integrate<F: Fn(&[f64]) -> f64>(&mut self, f: &F) -> IntegrateResult;
}
```

**功能**：数值积分器的统一 trait。被积函数 `f` 接收单位超立方体 $[0,1]^d$ 中的一个点，返回标量值。物理域上的积分需要在闭包内手动进行线性变量替换。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&F where F: Fn(&[f64]) -> f64` | 被积函数，输入为 $[0,1]^d$ 中的坐标点 |

**返回值**：`IntegrateResult`——积分估计值与标准误差。

**说明**：`integrate` 接收 `&F` 而非 `F`，因此可以多次调用同一闭包而无需重新分配。

---

## VegasOptions

**签名**：

```rust
#[derive(Debug, Clone, Copy)]
pub struct VegasOptions {
    pub n_bins: usize,
    pub n_samples: usize,
    pub iterations: usize,
    pub learning_rate: f64,
    pub seed: u64,
}
```

**功能**：Vegas 积分器的调优参数。

**字段**：

| 字段 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `n_bins` | `usize` | `64` | 每个维度的分箱数（总网格为各维 1-D 网格的直积） |
| `n_samples` | `usize` | `10_000` | 每次迭代的采样数 |
| `iterations` | `usize` | `10` | 自适应迭代次数 |
| `learning_rate` | `f64` | `1.5` | 网格平滑/学习率（典型范围 1.0–2.0；值越大，网格更新越平缓） |
| `seed` | `u64` | `0x0C45` | 随机数种子（跨运行确定性） |

**实现细节**：

- `Default` 实现提供上述默认值。
- 学习率的作用：在网格更新时，每个 bin 的重要性值 $d_i$ 被变换为 $d_i^{1/\text{learning\_rate}}$。`learning_rate = 1.0` 为完整步长，更大的值使网格变化更平滑、更稳定。
- 使用 `Xoshiro256PlusPlus` 伪随机数生成器（通过 `seed` 初始化），确保结果可复现。

**示例**：

```rust
use ocas_eval::numeric::VegasOptions;

// 高精度配置
let opts = VegasOptions {
    n_bins: 100,
    n_samples: 100_000,
    iterations: 20,
    learning_rate: 1.5,
    seed: 42,
};

// 使用默认值
let default_opts = VegasOptions::default();
```

**参见**：[`Vegas::new`](#vegasnew)

---

## Vegas

**签名**：

```rust
pub struct Vegas { /* 私有字段 */ }
```

**功能**：基于 Lepage Vegas 算法的自适应蒙特卡洛积分器，在单位超立方体 $[0,1]^d$ 上工作。

**设计原理**：Vegas 维护一个乘积网格（每个维度独立的 1-D 分箱），通过迭代精炼使每个 bin 捕获被积函数方差的近似相等份额。多轮迭代的估计值通过逆方差加权组合（见 [`StatisticsAccumulator`](#statisticsaccumulator)）。

### Vegas::new

**签名**：

```rust
pub fn new(n_dims: usize, opts: VegasOptions) -> Self
```

**功能**：创建一个 `n_dims` 维的 Vegas 积分器。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n_dims` | `usize` | 积分维度数（≥ 1） |
| `opts` | `VegasOptions` | 调优参数 |

**返回值**：初始化的 `Vegas` 实例，每个维度使用均匀分箱。

**示例**：

```rust
use ocas_eval::numeric::{Vegas, VegasOptions};

// 二维积分器，使用默认参数
let mut vegas = Vegas::new(2, VegasOptions::default());
```

---

### Vegas::integrate

**签名**：

```rust
pub fn integrate<F: Fn(&[f64]) -> f64>(&mut self, f: &F) -> IntegrateResult
```

**功能**：执行积分计算。对每个迭代：在单位超立方体上采样 `n_samples` 个点，计算被积函数值，累加统计量，然后更新网格边界以优化采样效率。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `&F where F: Fn(&[f64]) -> f64` | 被积函数，输入为 $[0,1]^d$ 中的坐标 |

**返回值**：`IntegrateResult`——所有迭代的逆方差加权平均估计及其标准误差。

**算法步骤**（每次迭代）：

1. **采样**：对每个维度，均匀选择一个 bin，再在 bin 内均匀采样；Jacobian = `n_bins × bin_width`（即逆 pdf 贡献）。
2. **求值**：计算 $f(\mathbf{x})$ 并以总 Jacobian 为权重累加到 [`StatisticsAccumulator`](#statisticsaccumulator)。
3. **训练信号**：将 $f^2 \cdot w$ 累加到对应 bin 的训练槽。
4. **统计收尾**：调用 `finalize_iteration()` 将本轮估计折叠到跨轮逆方差加权平均中。
5. **网格更新**：对每个维度，执行平滑累积弧长重分布——3-bin 平滑、均值归一化、阻尼 $d^{1/\text{lr}}$、在等弧长分位数重分布边界、单调性钳制。

**示例**：

```rust
use ocas_eval::numeric::{Vegas, VegasOptions};

// 计算 ∫₀¹∫₀¹ sin(x·y) dx dy
let opts = VegasOptions {
    n_samples: 50_000,
    iterations: 15,
    ..Default::default()
};
let mut vegas = Vegas::new(2, opts);
let r = vegas.integrate(&|x: &[f64]| (x[0] * x[1]).sin());
println!("结果：{:.6} ± {:.6}", r.integral, r.error);
// 结果：≈ 0.2398 ± 0.0002
```

**注意**：被积函数接收的是 $[0,1]^d$ 中的坐标。若需在物理域 $[a,b]$ 上积分，需在闭包内进行线性变换 $x_{\text{phys}} = a + u \cdot (b-a)$ 并乘以 Jacobian $(b-a)$。一维情况下可直接使用 [`integrate_1d`](#integrate_1d)。

---

### Vegas::result

**签名**：

```rust
pub fn result(&self) -> IntegrateResult
```

**功能**：获取最近一次 `integrate` 调用后的累积估计结果。

**返回值**：`IntegrateResult`——当前最佳积分估计和标准误差。

---

### Vegas::iterations

**签名**：

```rust
pub fn iterations(&self) -> usize
```

**功能**：返回已完成的迭代次数。

---

## integrate_1d

**签名**：

```rust
pub fn integrate_1d<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    opts: VegasOptions,
) -> IntegrateResult
```

**功能**：一维数值积分的便捷函数。在物理区间 $[a, b]$ 上积分函数 $f(x)$，内部自动处理变量替换 $u \mapsto a + u(b-a)$ 并将 Jacobian $(b-a)$ 折入结果。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `F where F: Fn(f64) -> f64` | 被积函数，直接接收物理坐标 $x \in [a, b]$ |
| `a` | `f64` | 积分下限 |
| `b` | `f64` | 积分上限（可小于 `a`，此时结果为负） |
| `opts` | `VegasOptions` | Vegas 参数（使用 `Default::default()` 获取默认值） |

**返回值**：`IntegrateResult`——积分估计值和标准误差。

**实现细节**：创建一维 `Vegas` 实例，将 `f` 包装为 `|u| f(a + u[0] * width) * width`，然后调用 `integrate`。

**示例**：

```rust
use ocas_eval::numeric::{integrate_1d, VegasOptions};

// ∫₀¹ x dx = 1/2
let r = integrate_1d(|x| x, 0.0, 1.0, Default::default());
assert!((r.integral - 0.5).abs() < 0.01);

// ∫₀² x² dx = 8/3 ≈ 2.667
let r = integrate_1d(|x| x * x, 0.0, 2.0, Default::default());
assert!((r.integral - 8.0 / 3.0).abs() < 0.05);

// 高精度配置
let opts = VegasOptions {
    n_samples: 100_000,
    iterations: 20,
    ..Default::default()
};
let r = integrate_1d(|x| (-x * x).exp(), -5.0, 5.0, opts);
// ∫₋₅⁵ e^{-x²} dx ≈ √π ≈ 1.7725
assert!((r.integral - std::f64::consts::PI.sqrt()).abs() < 0.01);
```

**参见**：[`Vegas::new`](#vegasnew)、[`Vegas::integrate`](#vegasintegrate)

---

## StatisticsAccumulator

**签名**：

```rust
#[derive(Debug, Clone)]
pub struct StatisticsAccumulator { /* 私有字段 */ }
```

**功能**：逆方差加权统计累加器，用于跨迭代合并蒙特卡洛积分估计。实现 Lepage 原始 Vegas 论文中的加权方案。

**设计原理**：每轮迭代产生一个积分估计 $I_i$ 和标准误差 $\sigma_i$。跨轮合并使用逆方差权重 $w_i = 1/\sigma_i^2$：

$$
I_{\text{combined}} = \frac{\sum_i w_i I_i}{\sum_i w_i}, \quad \sigma_{\text{combined}} = \frac{1}{\sqrt{\sum_i w_i}}
$$

同时维护 $\chi^2$ 统计量 $\sum_i w_i (I_i - I_{\text{combined}})^2$ 用于诊断网格分层质量——$\chi^2 \approx \text{iterations} - 1$ 表示各轮估计一致。

**内部字段**（私有）：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sum_w` | `f64` | 当前迭代：$\sum \text{weight}$ |
| `sum_wf` | `f64` | 当前迭代：$\sum \text{weight} \cdot f$ |
| `sum_wf2` | `f64` | 当前迭代：$\sum \text{weight} \cdot f^2$ |
| `integral` | `f64` | 跨迭代加权平均积分估计 |
| `error` | `f64` | 积分估计的标准误差 |
| `chi_square` | `f64` | 跨迭代 $\chi^2$ 统计量 |
| `iterations` | `usize` | 已完成迭代数 |

### StatisticsAccumulator::new

**签名**：

```rust
pub fn new() -> Self
```

**功能**：创建空累加器。初始状态：`integral = 0.0`，`error = INFINITY`，`iterations = 0`。

---

### StatisticsAccumulator::add_sample

**签名**：

```rust
pub fn add_sample(&mut self, weight: f64, f: f64)
```

**功能**：添加一个采样点。在 Vegas 中，`weight` 是该点的 Jacobian（逆 pdf），`f` 是被积函数值。该采样对本轮积分估计的贡献为 `weight * f`。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `weight` | `f64` | 采样权重（Vegas 中为 Jacobian = `n_bins × bin_width` 各维度之积） |
| `f` | `f64` | 被积函数在该点的值 |

**内部累加**：`sum_w += weight`，`sum_wf += weight * f`，`sum_wf2 += weight * f²`。

---

### StatisticsAccumulator::finalize_iteration

**签名**：

```rust
pub fn finalize_iteration(&mut self)
```

**功能**：收尾当前迭代——计算本轮均值和方差，通过逆方差加权折叠到跨轮平均中，然后重置本轮累加器。

**本轮均值**：$\bar{f} = \sum(w \cdot f) / \sum w$

**本轮方差**：$\text{Var} = \sum(w \cdot f^2) / \sum w - \bar{f}^2$

**退化处理**：若 `sum_w ≤ 0` 或方差为负（数值问题），跳过本轮但仍重置。

**误差钳制**：标准误差被钳制到不小于 $10^{-150}$（避免逆方差权重溢出）。

---

### StatisticsAccumulator::samples

**签名**：

```rust
pub fn samples(&self) -> usize
```

**功能**：返回当前（未收尾）迭代的近似采样数。

**注意**：内部通过 `sum_w` 推导（`sum_w as usize`），当权重为 1 时精确；Vegas 的 Jacobian 权重使此值为近似值。

---

### StatisticsAccumulator::integral

**签名**：

```rust
pub fn integral(&self) -> f64
```

**功能**：当前最佳积分估计（跨迭代逆方差加权平均）。

---

### StatisticsAccumulator::error

**签名**：

```rust
pub fn error(&self) -> f64
```

**功能**：积分估计的标准误差，即 $\sigma_{\text{combined}} = 1/\sqrt{\sum w_i}$。

---

### StatisticsAccumulator::chi_square

**签名**：

```rust
pub fn chi_square(&self) -> f64
```

**功能**：跨迭代的 $\chi^2$ 统计量。值接近 `iterations - 1` 表示各轮估计一致（网格分层良好）；显著偏大则表明被积函数存在未被网格捕获的结构。

---

### StatisticsAccumulator::iterations

**签名**：

```rust
pub fn iterations(&self) -> usize
```

**功能**：已完成的迭代次数。

---

**完整示例**：

```rust
use ocas_eval::numeric::StatisticsAccumulator;

let mut acc = StatisticsAccumulator::new();

// 模拟两轮迭代，每轮三个采样
// 第一轮
acc.add_sample(1.0, 1.5);
acc.add_sample(1.0, 1.8);
acc.add_sample(1.0, 2.1);
acc.finalize_iteration();
assert_eq!(acc.iterations(), 1);

// 第二轮
acc.add_sample(1.0, 1.6);
acc.add_sample(1.0, 1.9);
acc.add_sample(1.0, 2.0);
acc.finalize_iteration();
assert_eq!(acc.iterations(), 2);

// 跨轮加权平均
println!("积分：{:.4} ± {:.4}", acc.integral(), acc.error());
println!("χ² = {:.2}", acc.chi_square());
```

**参见**：[`Vegas::integrate`](#vegasintegrate)（内部使用此累加器）

---

## 与 ExpressionEvaluator 的配合

Vegas 积分器可与 [`ExpressionEvaluator`](./rust-evaluation.md) 配合，对符号表达式进行数值积分：

```rust
use ocas_eval::numeric::{Vegas, VegasOptions};
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("sin", &[ctx.var("x")]);

let evaluator: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();
let mut vegas = Vegas::new(1, VegasOptions::default());
let r = vegas.integrate(&|x: &[f64]| {
    evaluator.evaluate(&[x[0]]).unwrap()[0]
});
```

对于一维情况，也可直接使用 `integrate_1d` 包装：

```rust
use ocas_eval::numeric::integrate_1d;
use ocas_eval::ExpressionEvaluator;
use ocas_atom::AtomArena;
use ocas_core::arena::Arena;

let arena = Arena::new();
let ctx = AtomArena::new(&arena);
let expr = ctx.fun("sin", &[ctx.var("x")]);
let evaluator: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(expr).unwrap();

let r = integrate_1d(
    |x| evaluator.evaluate(&[x]).unwrap()[0],
    0.0,
    std::f64::consts::PI,
    Default::default(),
);
// ∫₀^π sin(x) dx ≈ 2.0
```

**参见**：[求值与 JIT](./rust-evaluation.md)、[蒙特卡洛积分数学基础](../math/monte-carlo-integration.md)
