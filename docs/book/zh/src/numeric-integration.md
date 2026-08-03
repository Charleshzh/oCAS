# 数值积分

oCAS 内置自适应蒙特卡洛积分器（Vegas），用于数值计算定积分。当符号反
导数不存在或计算代价过高时，数值积分是实用的替代方案。

---

## Vegas 算法

Vegas（Lepage 1978）是一种重要性采样蒙特卡洛方法，通过迭代细化分段常数
近似来集中采样到被积函数变化最剧烈的区域。oCAS 实现了标准自适应网格变
体，可配置迭代次数、采样预算和学习率。

| 入口 | 说明 |
|---|---|
| `integrate_1d(f, a, b, opts)` | 一维积分便捷函数 |
| `Vegas::new(n_dims, opts)` | 多维积分器 |
| `Integrator::integrate(&mut self, f)` | 运行积分 |

返回 `IntegrateResult { integral, error }`。

---

## 快速上手：一维积分

```rust
use ocas_eval::numeric::{integrate_1d, VegasOptions};

// 计算 sin(x) 从 0 到 π 的积分
let opts = VegasOptions {
    n_bins: 50,
    n_samples: 10_000,
    iterations: 5,
    learning_rate: 1.5,
    seed: 42,
};
let result = integrate_1d(|x| x.sin(), 0.0, std::f64::consts::PI, opts);
println!("integral = {:.6} ± {:.6}", result.integral, result.error);
// ≈ 2.000000
```

---

## 多维积分

对多变量积分，直接使用 `Vegas`：

```rust
use ocas_eval::numeric::{Vegas, VegasOptions, Integrator};

let opts = VegasOptions {
    n_bins: 50,
    n_samples: 10_000,
    iterations: 10,
    learning_rate: 1.5,
    seed: 0,
};
let mut vegas = Vegas::new(2, opts);

// 计算 x*y 在 [0,1]×[0,1] 上的积分 → 精确值 = 0.25
let result = vegas.integrate(&|coords| coords[0] * coords[1]);
println!("integral = {:.6} ± {:.6}", result.integral, result.error);
```

闭包接收一个 `&[f64]` 切片，每个元素对应一个维度，已映射到单位超立方
体 `[0, 1]ⁿ`。

---

## 调优参数

| 字段 | 默认值 | 作用 |
|---|---|---|
| `n_bins` | 64 | 每个维度的网格分箱数 |
| `n_samples` | 10,000 | 每次迭代的采样数 |
| `iterations` | 10 | 自适应迭代次数 |
| `learning_rate` | 1.5 | 网格自适应速度（1.0–2.0） |
| `seed` | `0x0C45`（`u64`） | 随机数种子，用于可重现性 |

增加迭代和采样次数可降低方差估计但增加运行时间。学习率控制网格对被积
函数结构的自适应强度。

---

## 统计累加器

`StatisticsAccumulator` 是 Vegas 内部使用的逆方差加权累加器。它跟踪每
次迭代的积分估计、卡方诊断和最终组合结果。

```rust
use ocas_eval::numeric::StatisticsAccumulator;

let mut acc = StatisticsAccumulator::new();
acc.add_sample(1.0, 1.5);
acc.add_sample(1.0, 1.8);
acc.add_sample(1.0, 2.1);
acc.finalize_iteration();

println!("integral = {:.6}", acc.integral());
println!("error    = {:.6}", acc.error());
println!("chi²     = {:.2}", acc.chi_square());
```

---

## Python 与 C 用法

### Python

```python
import ocas

# 一维便捷函数
result = ocas.integrate_1d(lambda x: x**2, 0, 1, n_samples=10000, iterations=10)
print(result.integral, result.error)

# 多维积分
vegas = ocas.Vegas(n_dims=2, n_samples=10000, iterations=10)
result = vegas.integrate(lambda coords: coords[0] * coords[1])
print(result.integral)
```

### C

```c
#include <ocas.h>

/* 一维便捷函数（opts 传 NULL 使用库默认参数） */
int err = 0;
struct ocas_OcasIntegrateResult result =
    ocas_integrate_1d(my_fn, NULL, 0.0, 1.0, NULL, &err);
printf("integral = %f ± %f\n", result.integral, result.error);
```

完整的绑定文档见 [Python API](./api/python.md) 和
[C/C++ API](./api/c.md) 章节。

---

## 限制

- Vegas 使用蒙特卡洛采样，结果是带误差条的统计估计，非精确符号值。
- 对高维且变量间强相关的积分，收敛较慢。
- 被积函数必须是纯 `f64 → f64` 函数；符号表达式需先编译为数值求值器
  （见 [求值与 JIT](./evaluation.md)）。
