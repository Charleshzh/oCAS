# 张量

oCAS 提供带显式指标管理的基础张量代数。`Tensor` 携带命名指标槽
（上标/下标）、可选的对称性元数据，并可转换为 `Atom` 表达式进行符号处理。

---

## 指标槽

每个张量指标具有：

| 属性 | 类型 | 含义 |
|---|---|---|
| `label` | `Atom` | 指标名（通常为希腊字母或 `i`、`j` 等） |
| `position` | `IndexPosition` | `Upper`（逆变）或 `Lower`（协变） |

```rust
use ocas_atom::tensor::{IndexSlot, IndexPosition};

let mu = IndexSlot::new("mu", IndexPosition::Upper);
let nu = IndexSlot::new("nu", IndexPosition::Lower);
```

---

## 创建张量

`Tensor` 由名称（`Symbol`）和指标槽向量定义：

```rust
use ocas_atom::tensor::{Tensor, IndexSlot, IndexPosition, Symmetry};

let slots = vec![
    IndexSlot::new("mu", IndexPosition::Upper),
    IndexSlot::new("nu", IndexPosition::Lower),
];
let t = Tensor::new(Symbol::new("g"), slots);

println!("rank = {}", t.rank());       // 2
println!("name = {}", t.name());       // g
```

可选地附加对称性元数据：

```rust
let symmetric = t.clone().with_symmetry(Symmetry::Symmetric);
let antisymmetric = t.clone().with_symmetry(Symmetry::Antisymmetric);
```

对称性元数据对下游消费者是建议性的；张量代数运算不会自动对称化。

---

## 缩并

`contract` 接受两个张量，返回标量（当所有指标都缩并时）或含剩余自由
指标的 `TensorProduct`：

```rust
use ocas_atom::tensor::{contract, Contracted};

// 缩并两个 1 阶张量：A^μ B_μ → 标量
let a = Tensor::new(Symbol::new("A"), vec![
    IndexSlot::new("mu", IndexPosition::Upper),
]);
let b = Tensor::new(Symbol::new("B"), vec![
    IndexSlot::new("mu", IndexPosition::Lower),
]);

match contract(&ctx, &a, &b) {
    Contracted::Scalar(expr) => println!("scalar: {}", expr),
    Contracted::Product(tp) => {
        for f in &tp.factors {
            println!("free factor: {}", f.name());
        }
    }
}
```

缩并按标签和相反位置匹配指标。单次调用可缩并多对匹配指标。

---

## 对称化

`symmetrise_sign` 返回张量反对称化的符号：对哑指标（已缩并标签）的
偶排列返回 `+1`，奇排列返回 `-1`。这在计算类行列式表达式时很有用。

```rust
use ocas_atom::tensor::symmetrise_sign;

let sign = symmetrise_sign(&tensor);
```

---

## 转换为 Atom

张量可通过 `to_atom` 降级为标准 `Atom` 表达式：

```rust
let atom = tensor.to_atom(&ctx);
println!("{}", atom);  // 例如 g(mu, nu)
```

这允许张量参与通用符号处理流水线（重写、化简、微分）。

---

## Python 与 C 用法

### Python

```python
from ocas import Tensor, contract_tensors, tensor_symmetrise_sign

# 创建 2 阶张量
g = Tensor("g", [("mu", "upper"), ("nu", "lower")])
print(g.rank())    # 2
print(g.name())    # g

# 缩并两个张量
A = Tensor("A", [("mu", "upper")])
B = Tensor("B", [("mu", "lower")])
result = contract_tensors(A, B)
print(result)      # 标量表达式

# 反对称化符号
sign = tensor_symmetrise_sign(g)
```

### C

```c
#include <ocas.h>

/* 创建张量 A^μ */
const char* labels[] = {"mu"};
int positions[] = {1};  /* 1 = upper */
ocas_OcasTensor* A = ocas_tensor_create("A", labels, positions, 1, &err);

/* 查询阶数 */
int rank = ocas_tensor_rank(A, &err);  /* 1 */

ocas_tensor_free(A);
```

完整的绑定文档见 [Python API](./bindings-python.md) 和
[C/C++ API](./bindings-c.md) 章节。

---

## 限制

- 张量代数是显式的：指标按标签匹配，不使用爱因斯坦求和约定。必须手
  动调用 `contract`。
- 基于图的规范化（自动对称性强制、Wick 缩并）推迟到 1.0 之后。
- 仅支持 `Symmetric` 和 `Antisymmetric` 对称类型，不支持混合或循环
  对称。
