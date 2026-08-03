# 高阶：张量代数与规范化

本章讨论 oCAS 中张量表达式的代数结构与规范化的数学基础。核心问题：给定一个包含缩并和对称性的张量乘积表达式，如何判定两个看起来不同的写法（例如 $T_{\mu\nu} U^{\nu\rho}$ 与 $U^{\rho\nu} T_{\nu\mu}$）是否代表同一个数学对象？oCAS 的方案是将张量表达式编码为有色有向图，然后通过图同构算法计算其**规范形式**。

---

## 前提知识

### 张量积与缩并

一个 $(r, s)$-型张量 $T$ 有 $r$ 个**逆变指标**（上标）和 $s$ 个**协变指标**（下标）。两个张量的**张量积** $T \otimes U$ 将指标简单拼接。

**缩并**（contraction）是将一个逆变指标和一个协变指标配对并求和的操作。在 Einstein 求和约定中：

$$
T_{\mu\nu} U^{\nu\rho} = \sum_{\nu} T_{\mu\nu} U^{\nu\rho}
$$

这里 $\nu$ 出现了两次——一次在下标、一次在上标——因此被求和掉，结果是一个 $(1,1)$-型张量。被求和掉的指标称为**哑指标**（dummy index），保留的指标称为**自由指标**（free index）。

> **Einstein 求和 vs. oCAS 显式匹配**：Einstein 约定隐式地对重复指标求和，依赖位置（上/下）判断配对。oCAS 采用显式匹配——每个指标槽是一个带位置标记的参数，缩并通过两个槽共享同一标签来表达。这消除了对"上/下位置"的隐式依赖，使系统更加明确。

### 指标升降

在有度规 $g_{\mu\nu}$ 的空间中，可以通过缩并度规来**升降**指标：

$$
T^{\mu} = g^{\mu\nu} T_{\nu}, \qquad T_{\mu} = g_{\mu\nu} T^{\nu}
$$

在 oCAS 的张量系统中，指标升降由 `IndexPosition` 枚举（`Upper` / `Lower`）标记在每个指标槽上，不自动引入度规张量。

---

## 基础概念

### 指标匹配

一个张量函数 $T(i, j, k)$ 的每个参数称为一个**指标槽**（index slot）。每个槽携带：

- **标签**（label）：一个 `Atom`，如 `i`、`mu`、`d0`
- **位置**（variance）：`Upper`（逆变）或 `Lower`（协变）

两个槽之间的**缩并**由标签匹配来表达：当两个槽共享同一标签且一个为 `Upper`、另一个为 `Lower` 时，它们形成一个缩并对。

### 哑指标管理

在一个张量乘积 $T(i, j) \cdot U(j, k)$ 中，标签 $j$ 出现恰好两次，因此是哑指标。oCAS 的哑指标管理遵循以下规则：

1. **识别**：遍历表达式中所有 `Fun` 节点的参数，统计每个标签的出现次数。出现恰好 **2** 次的标签是哑指标。
2. **验证**：规范化（`canon.rs`）中，标签出现超过 2 次报错 `TensorCanonError::ContractedMoreThanOnce`；`dummy.rs` 声明了 `OverContracted` 与 `BadContraction`（两次出现方差相同）错误变体，但 `refresh_dummies` 本身不执行方差校验。
3. **重命名**：规范化后，哑指标被替换为按组编号的新名称 `d0`, `d1`, `d2`, …。不同指标组（如时空指标和内部指标）使用不同的命名空间 `d{group}_{n}`。

```rust
// 哑指标重命名示例
// 输入：T(i, j) * U(j, k) * V(i, l)
// 输出：T(d0, d1) * U(d1, k) * V(d0, l)
// i → d0（组 0，第 0 个），j → d1（组 0，第 1 个）
```

这一重命名由 `dummy.rs` 中的 `refresh_dummies` 函数实现。它通过 `TensorRegistry::index_group` 查询标签所属的指标组，确保不同维度的哑指标使用不同的命名空间。

---

## 核心理论

### 图编码

张量规范化的第一步是将张量表达式编码为一个有色有向图 $G = (V, E)$。编码规则如下（实现于 `canon.rs` 的 `tensor_to_graph`）：

#### 顶点类型

| 类型 | 顶点数据 `TgNode` | 含义 |
|---|---|---|
| 头顶点 | `Head(hash)` | 一个张量函数头 $T$，`hash` 为函数名的哈希值 |
| 槽顶点 | `Slot(hash)` | 一个指标槽，`hash` 为指标标签的哈希值（对称槽统一为 `Slot(0)`） |
| 标量顶点 | `Scalar(hash)` | 乘积中的标量因子 |

#### 边类型

| 类型 | 边数据 `TgEdge` | 方向 | 含义 |
|---|---|---|---|
| 头→槽 | `HeadToSlot(pos, flag)` | 有向 | 张量头到其第 `pos` 个槽；`flag=1` 可见，`flag=0` 隐藏（对称） |
| 缩并 | `Contraction(group)` | 无向 | 两个缩并的槽顶点之间的边；`group` 为指标组标识 |

#### 对称性处理

当张量的某些槽是对称的（`SymmetrySpec::symmetric_subsets`）时，这些槽在图中被标记为**隐藏**：

- 对称槽的顶点数据统一为 `Slot(0)`，使得图同构引擎可以自由交换它们。
- 对称槽的头→槽边标记为 `HeadToSlot(pos, 0)`（flag=0），并使用**排序后的位置**而非原始位置，使 $T(a, b)$ 和 $T(b, a)$ 生成相同的编码（flag 本身仍参与证书比较）。
- 输入时，对称槽按标签字母序预排序，确保 $T(a, b)$ 和 $T(b, a)$ 生成相同的图。

**关键不变量**：两个张量表达式数学上等价 $\Leftrightarrow$ 它们的图同构。

### McKay 精炼-个体化算法

图同构判定和规范标号计算使用 McKay 精炼-个体化算法（nauty 系列），实现在 `graph.rs`（约 1100 行的独立实现）。

#### 算法概览

给定一个图 $G$，算法输出：

1. **规范标号** $\lambda: V \to \{0, 1, \dots, n-1\}$——一个双射，使得对任意同构图 $G'$，$\lambda$ 和 $\lambda'$ 产生的**证书**（certificate）完全相同。
2. **自同构轨道**——顶点在自同构群作用下的等价类。
3. **自同构群阶** $|\text{Aut}(G)|$。

#### 第 1 步：初始着色

将顶点按其 `data` 值分组：相同 `data` 的顶点属于同一**单元**（cell）。例如，所有 `Head(42)` 顶点在同一单元，所有 `Slot(7)` 在另一单元。

这定义了一个**有序划分** $\pi_0 = (C_1, C_2, \dots, C_k)$，其中每个 $C_i$ 是一个顶点集合。

#### 第 2 步：1-WL 颜色精炼

精炼过程迭代地将每个非平凡单元（$|C| > 1$）按其对其他单元的**邻居签名**（neighbour signature）进行分裂：

对于顶点 $v \in C_i$ 相对于目标单元 $C_j$ 的签名：

$$
\sigma(v, C_j) = \text{sorted}\bigl[\,(\text{edge\_data}, \text{direction}) \;\big|\; \text{edge } v \to w,\; w \in C_j\,\bigr]
$$

其中 `direction` 编码为：`0` = 有向出边，`1` = 有向入边，`2` = 无向边。

如果同一单元内的顶点有不同的签名，则将该单元按签名分裂为多个子单元。重复此过程直到划分是**公平的**（equitable）——即对于任意两个在同一单元的顶点 $u, v$ 和任意单元 $C_j$，$\sigma(u, C_j) = \sigma(v, C_j)$。

**复杂度**：每轮精炼最多产生 $O(n)$ 个新单元，每步需要 $O(n^2)$ 时间检查所有 (单元, 目标) 对。实际中，精炼通常在 $O(n \log n)$ 步内收敛。

#### 第 3 步：个体化-精炼搜索（DFS）

如果精炼后划分仍然非离散（存在 $|C| > 1$ 的单元），算法进入搜索阶段：

1. **选择最小非平凡单元** $C_{\min}$。
2. **个体化**（individualize）$C_{\min}$ 中的每个顶点 $v$：将 $v$ 从 $C_{\min}$ 移出，放入一个新的单元素单元（位于 $C_{\min}$ 之前）。这打破了 $v$ 与其他顶点的对称性。
3. 对个体化后的划分重新精炼。
4. 递归直到划分离散（每个单元恰好一个顶点），产生一个**离散标号**。

搜索使用**栈式 DFS**（非递归），通过 `SearchFrame` 记录每层的划分、路径不变量和是否在最左路径。

#### 第 4 步：路径不变量与自同构剪枝

##### 路径不变量

每个搜索节点记录一个**路径不变量** $\mathcal{I}$——从根到当前节点的划分单元长度序列：

$$
\mathcal{I} = \bigl[\,|C_1^{(0)}|, \dots, |C_k^{(0)}|\,\bigr] \;\to\; \bigl[\,|C_1^{(1)}|, \dots, |C_{k'}^{(1)}|\,\bigr] \;\to\; \cdots
$$

如果两个搜索节点的路径不变量不同，它们不可能产生同构的离散标号，因此可以安全剪枝。（当前实现计算并随 `SearchFrame` 记录该不变量，但尚未执行基于不变量的剪枝——实际生效的剪枝只有自同构轨道剪枝。）

##### 自同构轨道剪枝

当算法发现两个不同的离散标号产生相同的证书时，它们之间的映射构成一个自同构。这个自同构被记录为**轨道生成器**（orbit generator）。

在后续搜索中，如果待个体化的顶点 $v$ 与已处理的顶点 $w$ 在自同构群作用下处于同一轨道（即存在自同构将 $v$ 映射到 $w$），则跳过 $v$——因为个体化 $v$ 必然产生与个体化 $w$ 同构的子树。

轨道的计算通过 BFS 闭包实现：从种子顶点出发，反复应用生成器（及其逆元），直到闭合。

#### 第 5 步：规范形式

在搜索过程中，算法维护当前找到的**字典序最大证书**。每当到达一个离散节点，计算其证书：

$$
\text{cert}(\lambda) = \text{sorted}\bigl[\,(i, j, \text{edge\_data}, \text{dir}) \;\big|\; \text{edge } \lambda^{-1}(i) \to \lambda^{-1}(j)\,\bigr]
$$

如果新证书更大，则更新最佳标号；如果相等，则记录自同构生成器。

最终输出：

- **规范图**：按最佳标号重新标号的图。
- **自同构群阶**：通过 BFS 闭包枚举生成器产生的所有置换来计算（上限 $n!$，实际中远小于此）。
- **轨道划分**：由 `compute_orbits` 从生成器集合计算。

#### 完整示例

考虑表达式 $T_{\mu\nu} U^{\nu\rho}$，其中 $T$ 和 $U$ 无对称性：

**图编码**：

```
顶点: H_T (Head), S_mu (Slot), S_nu_T (Slot), H_U (Head), S_nu_U (Slot), S_rho (Slot)
边:
  H_T → S_mu    (HeadToSlot(0, 1))
  H_T → S_nu_T  (HeadToSlot(1, 1))
  H_U → S_nu_U  (HeadToSlot(0, 1))
  H_U → S_rho   (HeadToSlot(1, 1))
  S_nu_T — S_nu_U (Contraction(0))   // 无向
```

**初始着色**：顶点按 `data` 分组——`{H_T}` 与 `{H_U}` 是两个独立的 Head 单元（函数名哈希不同）；槽顶点按标签哈希分组：`{S_mu}`、`{S_nu_T, S_nu_U}`（标签同为 `nu`）、`{S_rho}`。

**精炼**：非平凡单元 `{S_nu_T, S_nu_U}` 对目标单元 `{H_T}` 的邻居签名不同（$S_{\nu_T}$ 有来自 $H_T$ 的入边——头→槽边从 $H_T$ 指向 $S_{\nu_T}$——而 $S_{\nu_U}$ 没有），因此被分裂。

**结果**：经过精炼后到达离散划分，无需搜索。规范形式确定了张量因子的顺序和哑指标的命名。

---

### Young 表与对称化投影

Young 表是表示对称群不可约表示的经典工具，在张量代数中用于构造具有特定对称性的张量分量。

#### Young 表的定义

一个 **Young 表**（Young tableau）由一个整数分拆 $\lambda = (\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_k > 0)$ 定义，其中 $\lambda_i$ 是第 $i$ 行的方格数。总方格数 $|\lambda| = \sum \lambda_i$ 等于张量的秩。

例如，$\lambda = (2, 1)$ 对应的 Young 图：

```
□ □
□
```

#### Young 对称化子

给定 Young 表 $\lambda$，其 **Young 对称化子** $e_\lambda$ 作用于张量 $T$ 产生具有 $\lambda$ 指定对称性的分量。定义如下：

1. **行对称化**：对每一行内的指标进行完全对称化。
2. **列反对称化**：对每一列内的指标进行完全反对称化。
3. **组合**：先列反对称化，再行对称化。

形式化地，设 $R$ 为行对称群的直积，$C$ 为列反对称群的直积，则：

$$
e_\lambda = \left(\sum_{\sigma \in R} \sigma\right) \cdot \left(\sum_{\tau \in C} \operatorname{sgn}(\tau) \cdot \tau\right)
$$

对张量 $T$ 的投影为：

$$
(e_\lambda \cdot T)_{i_1 \dots i_n} = \sum_{\sigma \in R} \sum_{\tau \in C} \operatorname{sgn}(\tau) \cdot T_{i_{\sigma(\tau(1))} \dots i_{\sigma(\tau(n))}}
$$

#### 置换的符号计算

对于给定的排列 $\pi$，Young 投影的符号计算分两步：

1. **行约束检查**：对每个位置 $p$，检查 $\pi(p)$ 是否与 $p$ 在同一行。如果不是（且该行长度 > 1），则 $\operatorname{sgn} = 0$（该排列不保持表的形状）。

2. **列奇偶性**：对每一列 $c$，提取 $\pi$ 限制在该列位置上的子置换，计算其奇偶性（通过圈分解：长度为偶数的圈贡献 $-1$ 因子）。所有列的奇偶性相乘得到总符号。

$$
\operatorname{sgn}_\lambda(\pi) = \begin{cases} \displaystyle\prod_{c} (-1)^{\#\text{even cycles in column } c} & \text{if } \pi \text{ preserves rows} \\ 0 & \text{otherwise} \end{cases}
$$

#### oCAS 实现：显式展开

`young.rs` 中的实现采用**显式置换展开**（非 BSGS 群论方法）：

1. 将行组 $R$（每行内部置换）与列组 $C$（每列内部置换）分别展开为完整位置映射，共 $\prod r_i! \cdot \prod c_j!$ 个组合。
2. 对每个组合 $(\sigma, \tau)$ 构造复合置换 $\tau \circ \sigma$，项系数为列置换奇偶 $\operatorname{sgn}(\tau)$。
3. 每项贡献 $\operatorname{sgn}(\tau) \cdot T(i_{\tau\sigma(1)}, \dots, i_{\tau\sigma(n)})$。
4. 所有项求和——这是经典 Young 对称子 $c_\lambda = a_\lambda b_\lambda$，**不做** Hook 长度归一化，每项系数仅为 $\pm 1$；因此结果是模标量 $\kappa = (\prod r_i! \cdot \prod c_j!) / \dim(\lambda)$ 的投影算子（$c_\lambda^2 = \kappa \cdot c_\lambda$）。

**完全反对称表** $\lambda = (1, 1, \dots, 1)$：行组平凡，列组为整个 $S_n$，产生标准交错和。

**完全对称表** $\lambda = (n)$：每列一个元素，列反对称化退化为恒等，行对称化产生所有排列的和。

---

### 对称性规格

`SymmetrySpec` 结构体声明一个张量函数头的槽位对称性，是图编码和 Young 投影之间的桥梁。

#### 三种对称类型

| 类型 | 字段 | 图编码行为 | Young 投影行为 |
|---|---|---|---|
| 对称子集 | `symmetric_subsets` | 槽标记为隐藏（`is_slot_hidden = true`），图同构可自由置换 | 保持对应行对称 |
| 反对称子集 | `antisymmetric_subsets` | 槽可见，编码与普通槽相同（无特殊标记；反对称符号由 Young 投影在表达式层面处理） | 保持对应列反对称 |
| 循环置换 | `cyclic` | 槽标记为隐藏 | 循环置换的对称化 |

#### `is_slot_hidden` 判定

一个槽位 $p$ 被标记为隐藏当且仅当：

$$
\text{is\_hidden}(p) \iff \exists S \in \text{symmetric\_subsets}: p \in S \;\lor\; \exists C = \text{cyclic}: p \in C
$$

隐藏槽在图中使用统一的颜色 `Slot(0)`，使图同构引擎可以自由交换它们，从而保证 $T(a, b) = T(b, a)$ 对对称张量成立。

#### TensorRegistry

`TensorRegistry` 管理所有张量的规格和指标组：

- `register(name, spec)`：注册张量函数头及其对称性规格。
- `set_index_group(label, group)`：设置指标标签所属的组（如时空指标组 1，内部指标组 2）。
- `index_group(label)`：查询指标组（0 为默认/未分组）。

指标组的作用是在哑指标重命名时防止跨维度的冲突。例如，时空哑指标 $\mu$ 和内部哑指标 $i$ 使用不同的命名空间 `d1_0`, `d1_1`, … 和 `d2_0`, `d2_1`, …。

---

## 在 oCAS 中的实现

### 模块结构

```
ocas-atom/src/tensor/
├── canon.rs      # 图编码 + 规范化入口
├── graph.rs      # McKay 精炼-个体化引擎（约 1100 行）
├── young.rs      # Young 表投影
├── spec.rs       # SymmetrySpec + TensorRegistry
└── dummy.rs      # 哑指标管理
```

### `canon.rs`：图编码与规范化

**入口函数** `canonicalize_tensors(ctx, expr, registry)`：

1. 如果表达式是加法（`Add`），逐项规范化并验证自由指标一致性。
2. 单项式进入 `canonicalize_single_term`：
   - **快速路径**：如果张量有对称子集、无反对称子集、且全部槽都是隐藏的（`is_slot_hidden` 全部为真），直接按标签字母序排序参数，无需图同构。
   - **通用路径**：调用 `tensor_to_graph` 编码为图，调用 `g.canonize()` 计算规范形式，调用 `reconstruct` 重建表达式。

**`tensor_to_graph`** 的编码逻辑：

- 对乘积中的每个因子 `Fun(name, args)`：
  - 创建头顶点 `Head(hash(name))`。
  - 对称槽按标签字母序预排序，使 $T(a,b)$ 和 $T(b,a)$ 生成相同的图。
  - 对每个参数创建槽顶点，隐藏槽使用 `Slot(0)`，可见槽使用 `Slot(hash(label))`。
  - 添加有向边 `Head → Slot`。
- 对乘积中标签出现恰好 2 次的槽对，添加无向缩并边。

**`reconstruct`** 的重建逻辑：

- 遍历规范图中的头顶点，按规范顺序收集张量因子。
- 对每个头的出边（按 `orig_pos` 或标签序排序），确定槽的标签：
  - 缩并槽：分配规范哑指标名 `d0`, `d1`, …
  - 自由槽：保留原始标签。
- 将因子组装为乘积或加法表达式。

### `graph.rs`：McKay 算法

这是 oCAS 中最大、最独立的算法模块（约 1100 行），实现了完整的 nauty 风格图同构引擎。

**数据结构**：

- `Graph<N, H, E>`：三参数泛型图。`N` = 顶点可见数据（参与比较），`H` = 顶点隐藏数据（不参与比较），`E` = 边数据（参与比较）。
- `Partition`：有序划分，`cells[i]` 为第 $i$ 个单元的顶点列表，`cell_of[v]` 为顶点 $v$ 所在单元的索引。
- `SearchFrame`：DFS 搜索栈帧，包含划分、路径不变量和是否在最左路径。
- `CanonicalForm<N, H, E>`：输出结果，包含 `vertex_map`（原→规范映射）、`orbits`（自同构轨道）、`automorphism_group_size`（群阶）和 `graph`（规范图）。

**关键实现细节**：

- **精炼**（`Partition::refine`）：使用 `stable_below` 指针优化——已稳定的单元不再参与分裂检查。当某单元被分裂后，`stable_below` 回退到该单元的索引，确保后续精炼覆盖所有受影响的单元。
- **邻居签名**（`cell_signatures` / `vertex_signature`）：对每个顶点 $v$，收集它到目标单元 $C_j$ 中所有邻居的 `(edge_data, direction)` 对，排序后作为签名。`direction` 区分有向出边（0）、有向入边（1）和无向边（2）。
- **证书**（`Graph::certificate`）：将图的所有边按规范标号转换为 `(i, j, data, dir)` 四元组列表，排序后作为字典序比较的键。
- **自同构群阶**（`group_size`）：通过 BFS 闭包枚举生成器产生的所有置换（**全枚举**，无截断；上限为 $n!$，实际远小于此）。
- **轨道计算**（`compute_orbits`）：BFS 从每个未访问顶点出发，沿生成器及其逆元扩展，直到闭合。

### `young.rs`：Young 投影

**`YoungTableau`** 结构体：

- `row_lengths: Vec<usize>`：各行的方格数，如 `[2, 1]` 表示 $\lambda = (2,1)$。
- `total_boxes()`：总方格数。
- `sign_of_permutation(perm)`：计算排列在给定表下的符号（$-1, 0, +1$）。

**`young_project(ctx, tensor_expr, tableau)`**：

- 仅处理 `Fun(name, args)` 形式的表达式。
- 使用 Heap 算法生成全排列（$O(n!)$ 时间，常数空间）。
- 对每个排列计算符号，非零则构造带符号的项。
- 返回所有项的加法（或单一项、零）。

**注意**：这是显式展开，适用于低秩张量（$n \leq 8$ 左右）。对于高秩张量，应考虑使用 BSGS（基本 Schreier–Sims）方法，但 oCAS 当前未实现。

**已知限制**：`sign_of_permutation` 的列奇偶性代码对**保持行但跨列移动**的置换会触发 `unwrap()` panic（在另一列的列位置查找失败）。完全对称表 $(n)$ 与完全反对称表 $(1, \ldots, 1)$ 不会触发；一般形状（如 $(2,1)$）的显式展开可能 panic，使用时需注意。

### `spec.rs`：对称性规格

**`SymmetrySpec`** 的三种构造器：

- `none()`：无对称性。
- `fully_symmetric(rank)`：所有槽完全对称。
- `fully_antisymmetric(rank)`：所有槽完全反对称。

**`is_slot_hidden(pos)`**：检查槽是否属于对称子集或循环子集。隐藏槽不参与图规范化的比较，允许图同构引擎自由置换它们。

### `dummy.rs`：哑指标管理

**`refresh_dummies(ctx, expr, registry)`**：

1. 收集表达式中所有 `Fun` 节点参数的出现次数。
2. 识别哑指标（出现恰好 2 次）。
3. 对每个哑指标，查询其指标组（`registry.index_group`），分配组内编号。
4. 重命名：组 0 的哑指标命名为 `d0`, `d1`, …；组 $g$ 的命名为 `d{g}_0`, `d{g}_1`, …。

**错误变体**（`dummy.rs` 声明，但 `refresh_dummies` 当前不触发）：

- `OverContracted(Symbol)`：指标出现超过 2 次。
- `BadContraction(Symbol)`：两次出现的方差相同（非上/下配对）。

（规范化路径中，超过 2 次的标签由 `canon.rs` 报 `TensorCanonError::ContractedMoreThanOnce`。）

---

## 进阶话题

### 隐藏槽与快速路径

对于全对称张量（如度规 $g_{\mu\nu}$），所有槽都被标记为隐藏。此时图中所有槽顶点共享同一颜色 `Slot(0)`，McKay 算法的自同构群将非常大（$|S_n| = n!$）。

为避免这种开销，`canonicalize_single_term` 检测到全对称张量时走**快速路径**：直接按标签字母序排序参数，无需构建图。这将规范化复杂度从 $O(n!)$ 降低到 $O(n \log n)$。

### 自同构群信息的用途

McKay 算法不仅输出规范形式，还输出自同构群的轨道和阶。这些信息在张量代数中有重要应用：

- **轨道**：标识哪些槽在对称性下等价，可用于化简求和。
- **群阶**：$|\text{Aut}(G)|$ 是表达式内在对称性的度量，可用于估计对称化展开的规模。

### 与 nauty 的关系

`graph.rs` 是一个独立的 nauty 风格实现，不依赖外部库。它实现了 McKay (1981/2014) 算法的核心组件：

- 1-WL 颜色精炼（等价于 nauty 的 `refine`）
- 个体化-精炼 DFS（等价于 `search`）
- 路径不变量记录（等价于 nauty 的 `compare`，但基于它的剪枝当前未启用）
- 自同构轨道剪枝（等价于 nauty 的 `orbits`）

与原始 nauty 的主要差异：

- 使用邻接表而非位矩阵，适合稀疏图。
- 证书使用排序的边列表而非规范邻接矩阵。
- 隐藏数据（`H` 参数）允许存储不参与比较的元数据。

### 扩展到高阶对称性

对于更复杂的对称性（如黎曼曲率张量 $R_{abcd}$ 的代数 Bianchi 恒等式），当前系统通过组合 `symmetric_subsets` 和 `antisymmetric_subsets` 来近似。完整的 Young 对称化通过 `young_project` 函数实现，但它产生的是展开后的和，而非保持压缩表示。

---

## 参考文献

1. **Cvitanović, P.** *Group Theory: Birdtracks, Lie's, and Exceptional Groups.* Princeton University Press, 2008. — 张量不变量与对称群表示的经典参考。
2. **McKay, B. D.** *Practical Graph Isomorphism.* Congressus Numerantium, 30:45–87, 1981 (updated 2014). — nauty 算法的原始论文。
3. **Faugère, J.-C.** *A new efficient algorithm for computing Gröbner bases (F4).* Journal of Pure and Applied Algebra, 139(1–3):61–88, 1999. — F4 算法，与张量缩并的矩阵方法相关。
4. **Cox, D., Little, J., O'Shea, D.** *Ideals, Varieties, and Algorithms.* Springer, 4th ed., 2015. — 多项式代数基础，Gröbner 基理论参考。
5. **Hamermesh, M.** *Group Theory and Its Application to Physical Problems.* Dover, 1989. — Young 表与张量对称性的物理应用。
