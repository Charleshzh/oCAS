# 自动微分

oCAS 通过*超对偶数*提供前向自动微分。这与符号微分（`diff`）不同：自
动微分在求函数值的同时数值计算导数，无需构建或化简表达式树。

---

## 超对偶数

超对偶数将标量域扩展出无穷小分量 ε₁, ε₂, …，满足 εᵢ² = 0。对这些
数的算术运算会精确传播一阶（以及可选的更高阶）偏导数。

oCAS 用 `HyperDual<T>` 表示超对偶数，其中 `T` 是实现 `DualCoeff` trait
的系数类型。实践中 `T` 通常是 `Rational`（标准构建中唯一完整实现
`DualCoeff` 的类型）。

| 类型 | 作用 |
|---|---|
| `DualShape` | 布局描述：变量分组、乘法表 |
| `HyperDual<T>` | 具体的对偶数，含值和导数槽 |
| `DualCoeff` | trait：`zero`、`one`、`T` 的算术运算 |
| `new_first_order(n)` | 便捷函数：含 `n` 个一阶分量的形状 |

---

## 快速上手

```rust
use ocas_domain::Rational;
use ocas_domain::dual::{new_first_order, HyperDual};

// 2 个变量的一阶对偶数
let shape = new_first_order::<Rational>(2);

// x = 1 + ε₁（变量 0，单位系数）
let x = HyperDual::variable(&shape, 0, Rational::new(1, 1));
// y = 2 + ε₂（变量 1，单位系数）
let y = HyperDual::variable(&shape, 1, Rational::new(2, 1));

// f = x * y → 值 = 2，∂f/∂x = 2，∂f/∂y = 1
let f = x * y;
println!("f    = {}", f.value());           // 2
println!("df/dx = {}", f.deriv(0).unwrap()); // 2
println!("df/dy = {}", f.deriv(1).unwrap()); // 1
```

---

## 支持的运算

`HyperDual<T>` 支持标准算术 trait（`Add`、`Sub`、`Mul`、`Div`、`Neg`）。
每个运算通过 `DualShape` 中预计算的乘法表传播导数。

| 运算 | 支持？ |
|---|---|
| `+`、`-`、`*`、`/` | 是 |
| `inv()`（乘法逆元） | 是 |
| `pow` / `exp` / `log` / 三角函数 | **否**（超越函数需要尚未实现的 `DualCoeff` 扩展） |

对超越函数，可回退到符号微分（`diff`）或将表达式编译为带自动微分支
持的数值求值器（计划在 1.0 之后实现）。

---

## 高阶导数

`DualShape` 支持多指标分组以计算高阶导数（如二阶 εᵢεⱼ 分量）。可使
用 `DualShape::new` 传入嵌套索引向量来自定义形状，或使用一阶便捷函数
处理常见场景。

---

## Python 与 C 用法

### Python

```python
from ocas import DualShape, HyperDual

shape = DualShape.first_order(2)
x = HyperDual.variable(shape, 0, 1)
y = HyperDual.variable(shape, 1, 2)

f = x * y
print(f.value())    # 2
print(f.deriv(0))   # 2  (∂f/∂x)
print(f.deriv(1))   # 1  (∂f/∂y)

# dunder 方法：__add__、__sub__、__mul__、__truediv__、__neg__
g = x + y * x
print(g.deriv(0))   # 3  (1 + y = 1 + 2)
```

### C

```c
#include <ocas.h>

int err = 0;
struct ocas_OcasDualShape *shape = ocas_dual_shape_new(2, &err);
struct ocas_OcasHyperDual *x = ocas_dual_variable(shape, 0, "1", &err);
struct ocas_OcasHyperDual *y = ocas_dual_variable(shape, 1, "2", &err);

struct ocas_OcasHyperDual *f = ocas_dual_mul(x, y, &err);
char *val = ocas_dual_value(f, &err);    /* "2" */
char *dx  = ocas_dual_deriv(f, 0, &err); /* "2" */

ocas_string_free(val);
ocas_string_free(dx);
ocas_hyperdual_free(f);
ocas_hyperdual_free(y);
ocas_hyperdual_free(x);
ocas_dual_shape_free(shape);
```

完整的绑定文档见 [Python API](./api/python.md) 和
[C/C++ API](./api/c.md) 章节。

---

## 限制

- 标准构建仅支持 `Rational` 系数（尚不支持浮点或有限域对偶数）。
- 超越函数（`sin`、`cos`、`exp`、`log`、`pow`）尚不可用于 `HyperDual`。
- 对偶表达式的 JIT 编译尚未集成。
