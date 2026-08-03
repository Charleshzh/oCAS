# 矩阵

> 来源：`ocas-poly/src/matrix.rs`

`Matrix<D>` 是 oCAS 中通用的稠密矩阵类型，元素来自任意 `EuclideanDomain`。采用行优先存储，所有消元操作使用 **无分数（fraction-free）** 算法，避免域元素的分数膨胀。

---

## `MatrixError`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    ShapeMismatch,
    RightHandSideIsNotVector,
    Inconsistent,
    Underdetermined { rank: usize },
    ResultNotInDomain,
}
```

矩阵运算中可能产生的错误。实现了 `Display` 和 `std::error::Error`。

### `ShapeMismatch`

**功能**：矩阵形状不兼容。

**触发场景**：`trace` 用于非方阵；`matmul` 的 `self.ncols != other.nrows`；`solve` 中 `self.nrows != b.len()`；`augment` 中行数不一致。

### `RightHandSideIsNotVector`

**功能**：右端不是列向量。

**触发场景**：当前版本中未直接由 `solve` 使用（`solve` 接受 `&[D::Element]` 切片），但保留供扩展。

### `Inconsistent`

**功能**：线性方程组无解（不相容）。

**触发场景**：`solve` 中消元后增广矩阵的自由行在最后一列（常数项）非零。

### `Underdetermined { rank }`

**功能**：线性方程组有无穷多解（欠定）。

**触发场景**：`solve` 中系数矩阵的秩 `rank` 小于变量数 `nvars`。

### `ResultNotInDomain`

**功能**：解不在目标域中。

**触发场景**：`solve` 的最终除法步骤中，`div` 返回 `None`（例如整数域中不能整除）。

**Display 输出**：

| 变体 | 输出 |
|---|---|
| `ShapeMismatch` | `"matrix shape mismatch"` |
| `RightHandSideIsNotVector` | `"right-hand side must be a column vector"` |
| `Inconsistent` | `"inconsistent linear system"` |
| `Underdetermined { rank }` | `"underdetermined system (rank N)"` |
| `ResultNotInDomain` | `"solution does not lie in the expected domain"` |

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;

// 不相容系统：两行相同但常数项不同
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(1)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(1), Integer::from(2)];
assert_eq!(a.solve(&b), Err(MatrixError::Inconsistent));
```

**参见**：[solve](#solve)

---

## `Matrix<D>`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<D: EuclideanDomain> {
    data: Vec<D::Element>,   // 行优先存储
    nrows: usize,
    ncols: usize,
    domain: D,
}
```

**功能**：域 `D` 上的稠密矩阵，元素以行优先顺序存储于 `Vec<D::Element>` 中。

**设计不变量**：
- `data.len() == nrows * ncols`（由构造函数保证）
- 元素 `(i, j)` 在 `data` 中的位置为 `i * ncols + j`
- `D` 约束为 `EuclideanDomain`，支持 `div`、`gcd` 等操作，使无分数消元成为可能

**Trait 实现**：
- `Index<(usize, usize)>` / `IndexMut<(usize, usize)>` — 通过 `matrix[(i, j)]` 访问元素
- `Display` — 每行一行，元素间空格分隔
- `Debug`、`Clone`、`PartialEq`、`Eq`

---

## 构造方法

### `new`

**签名**：`pub fn new(nrows: usize, ncols: usize, data: Vec<D::Element>, domain: D) -> Self`

**功能**：从行优先的一维数据创建矩阵。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `nrows` | `usize` | 行数 |
| `ncols` | `usize` | 列数 |
| `data` | `Vec<D::Element>` | 行优先排列的元素，长度必须等于 `nrows * ncols` |
| `domain` | `D` | 系数域 |

**返回值**：`Matrix<D>`

**错误**：若 `data.len() != nrows * ncols`，运行时 panic。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// 创建 2×3 矩阵：
// [ 1  2  3 ]
// [ 4  5  6 ]
let m = Matrix::new(2, 3, vec![
    Integer::from(1), Integer::from(2), Integer::from(3),
    Integer::from(4), Integer::from(5), Integer::from(6),
], d);
assert_eq!(m.nrows(), 2);
assert_eq!(m.ncols(), 3);
assert_eq!(m[(0, 0)], Integer::from(1));
assert_eq!(m[(1, 2)], Integer::from(6));
```

**参见**：[from_rows](#from_rows)、[zeros](#zeros)、[identity](#identity)

---

### `zeros`

**签名**：`pub fn zeros(nrows: usize, ncols: usize, domain: D) -> Self`

**功能**：创建指定形状的零矩阵。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `nrows` | `usize` | 行数 |
| `ncols` | `usize` | 列数 |
| `domain` | `D` | 系数域 |

**返回值**：所有元素为 `domain.zero()` 的 `Matrix<D>`。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let z = Matrix::zeros(2, 3, d);
assert_eq!(z[(0, 0)], Integer::from(0));
assert_eq!(z[(1, 2)], Integer::from(0));
```

**参见**：[identity](#identity)、[new](#new)

---

### `identity`

**签名**：`pub fn identity(n: usize, domain: D) -> Self`

**功能**：创建 `n × n` 单位矩阵。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `usize` | 矩阵维度 |
| `domain` | `D` | 系数域 |

**返回值**：对角线为 `domain.one()`、其余为 `domain.zero()` 的方阵。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let id = Matrix::identity(3, d);
assert_eq!(id[(0, 0)], Integer::from(1));
assert_eq!(id[(0, 1)], Integer::from(0));
assert_eq!(id[(2, 2)], Integer::from(1));
```

**参见**：[zeros](#zeros)

---

### `from_rows`

**签名**：`pub fn from_rows(rows: Vec<Vec<D::Element>>, domain: D) -> Self`

**功能**：从嵌套行向量创建矩阵。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `rows` | `Vec<Vec<D::Element>>` | 每个内层 `Vec` 是一行，所有行长度必须相同 |
| `domain` | `D` | 系数域 |

**返回值**：`Matrix<D>`

**错误**：若行长度不一致，运行时 panic。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![
        vec![Integer::from(1), Integer::from(2), Integer::from(3)],
        vec![Integer::from(4), Integer::from(5), Integer::from(6)],
    ],
    d,
);
assert_eq!(m.nrows(), 2);
assert_eq!(m.ncols(), 3);
assert_eq!(m[(0, 2)], Integer::from(3));
```

**参见**：[new](#new)、[into_rows](#into_rows)

---

## 访问器方法

### `nrows`

**签名**：`pub fn nrows(&self) -> usize`

**功能**：返回矩阵行数。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::new(3, 2, vec![Integer::from(0); 6], d);
assert_eq!(m.nrows(), 3);
```

---

### `ncols`

**签名**：`pub fn ncols(&self) -> usize`

**功能**：返回矩阵列数。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::new(3, 2, vec![Integer::from(0); 6], d);
assert_eq!(m.ncols(), 2);
```

---

### `shape`

矩阵没有单独的 `shape()` 方法。形状通过组合 `nrows()` 和 `ncols()` 获取：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::new(3, 4, vec![Integer::from(0); 12], d);
let (rows, cols) = (m.nrows(), m.ncols());
assert_eq!((rows, cols), (3, 4));
```

---

### `domain`

**签名**：`pub fn domain(&self) -> &D`

**功能**：返回矩阵系数域的引用。

**示例**：

```rust
use ocas_domain::{IntegerDomain, EuclideanDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::zeros(1, 1, d);
let dom = m.domain();
assert!(dom.is_zero(&dom.zero()));
```

---

### `data`

**签名**：`pub fn data(&self) -> &[D::Element]`

**功能**：返回底层行优先存储的只读切片。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let raw = m.data();
assert_eq!(raw.len(), 4);
assert_eq!(raw[0], Integer::from(1));  // (0,0)
assert_eq!(raw[3], Integer::from(4));  // (1,1)
```

---

### `row`

**签名**：`pub fn row(&self, i: usize) -> Vec<D::Element>`

**功能**：返回第 `i` 行的拷贝。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `i` | `usize` | 行索引，从 0 开始 |

**返回值**：`Vec<D::Element>`，长度等于 `ncols`。

**错误**：索引越界时 panic。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let r = m.row(1);
assert_eq!(r, vec![Integer::from(3), Integer::from(4)]);
```

**参见**：[column](#column)、[into_rows](#into_rows)

---

### `column`

**签名**：`pub fn column(&self, j: usize) -> Vec<D::Element>`

**功能**：返回第 `j` 列的拷贝。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `j` | `usize` | 列索引，从 0 开始 |

**返回值**：`Vec<D::Element>`，长度等于 `nrows`。

**错误**：索引越界时 panic。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let c = m.column(0);
assert_eq!(c, vec![Integer::from(1), Integer::from(3)]);
```

**参见**：[row](#row)

---

### `into_rows`

**签名**：`pub fn into_rows(self) -> Vec<Vec<D::Element>>`

**功能**：将矩阵消费并转换为嵌套行向量。

**返回值**：外层 `Vec` 长度为 `nrows`，每个内层 `Vec` 长度为 `ncols`。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let rows = m.into_rows();
assert_eq!(rows.len(), 2);
assert_eq!(rows[0], vec![Integer::from(1), Integer::from(2)]);
assert_eq!(rows[1], vec![Integer::from(3), Integer::from(4)]);
```

**参见**：[from_rows](#from_rows)、[row](#row)

---

## 索引与显示

### `Index<(usize, usize)>`

**签名**：`impl Index<(usize, usize)> for Matrix<D>`

**功能**：通过 `matrix[(i, j)]` 只读访问元素 `(i, j)`，等价于 `data[i * ncols + j]`。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(5), Integer::from(6)]], d,
);
assert_eq!(m[(0, 0)], Integer::from(5));
assert_eq!(m[(0, 1)], Integer::from(6));
```

---

### `IndexMut<(usize, usize)>`

**签名**：`impl IndexMut<(usize, usize)> for Matrix<D>`

**功能**：通过 `matrix[(i, j)] = val` 可写访问元素。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let mut m = Matrix::zeros(2, 2, d);
m[(0, 1)] = Integer::from(42);
m[(1, 0)] = Integer::from(-7);
assert_eq!(m[(0, 1)], Integer::from(42));
assert_eq!(m[(1, 0)], Integer::from(-7));
```

---

### `Display`

**签名**：`impl Display for Matrix<D> where D::Element: Display`

**功能**：格式化输出矩阵，每行占一行，元素间空格分隔。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let s = format!("{:?}", m);
// 输出：Debug 表示（IntegerDomain 未实现 Display，Matrix 的 Display 需要元素实现 Display）
assert!(s.contains("1"));
```

---

## 矩阵运算

### `swap_rows`

**签名**：`pub fn swap_rows(&mut self, i: usize, j: usize, start_col: usize)`

**功能**：交换第 `i` 行和第 `j` 行，从第 `start_col` 列开始。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `i` | `usize` | 第一行索引 |
| `j` | `usize` | 第二行索引 |
| `start_col` | `usize` | 起始列（含），此列之前的元素不变 |

**返回值**：无（就地修改）。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let mut m = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
m.swap_rows(0, 1, 0);
assert_eq!(m[(0, 0)], Integer::from(3));
assert_eq!(m[(0, 1)], Integer::from(4));
assert_eq!(m[(1, 0)], Integer::from(1));
assert_eq!(m[(1, 1)], Integer::from(2));
```

---

### `transpose`

**签名**：`pub fn transpose(&self) -> Matrix<D>`

**功能**：返回矩阵的转置。$(A^T)_{ij} = A_{ji}$。

**返回值**：新的 `Matrix<D>`，形状为 `(ncols, nrows)`。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2), Integer::from(3)],
         vec![Integer::from(4), Integer::from(5), Integer::from(6)]],
    d,
);
let t = a.transpose();
assert_eq!(t.nrows(), 3);
assert_eq!(t.ncols(), 2);
assert_eq!(t[(0, 0)], Integer::from(1));
assert_eq!(t[(0, 1)], Integer::from(4));
assert_eq!(t[(2, 1)], Integer::from(6));
```

**参见**：[trace](#trace)

---

### `trace`

**签名**：`pub fn trace(&self) -> Result<D::Element, MatrixError>`

**功能**：计算方阵的迹，即对角线元素之和：$\text{tr}(A) = \sum_{i=0}^{n-1} a_{ii}$。

**返回值**：
- `Ok(sum)` — 迹的值
- `Err(MatrixError::ShapeMismatch)` — 矩阵非方阵

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
// tr = 1 + 4 = 5
assert_eq!(a.trace().unwrap(), Integer::from(5));
```

**非方阵错误**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let m = Matrix::from_rows(vec![vec![Integer::from(1), Integer::from(2)]], d);
assert_eq!(m.trace(), Err(MatrixError::ShapeMismatch));
```

**参见**：[determinant](#determinant)

---

### `matmul`

**签名**：`pub fn matmul(&self, other: &Matrix<D>) -> Result<Matrix<D>, MatrixError>`

**功能**：计算矩阵乘积 $C = AB$，其中 $C_{ij} = \sum_k A_{ik} \cdot B_{kj}$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `other` | `&Matrix<D>` | 右乘矩阵 |

**返回值**：
- `Ok(Matrix)` — 乘积矩阵，形状为 `(self.nrows, other.ncols)`
- `Err(MatrixError::ShapeMismatch)` — `self.ncols != other.nrows`

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
let b = Matrix::from_rows(
    vec![vec![Integer::from(5), Integer::from(6)],
         vec![Integer::from(7), Integer::from(8)]],
    d,
);
let c = a.matmul(&b).unwrap();
// C = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]]
//   = [[19, 22], [43, 50]]
assert_eq!(c[(0, 0)], Integer::from(19));
assert_eq!(c[(0, 1)], Integer::from(22));
assert_eq!(c[(1, 0)], Integer::from(43));
assert_eq!(c[(1, 1)], Integer::from(50));
```

**形状不兼容**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let a = Matrix::from_rows(vec![vec![Integer::from(1), Integer::from(2)]], d);
let b = Matrix::from_rows(vec![vec![Integer::from(3), Integer::from(4)]], d);
// a: 1×2, b: 1×2 → ncols(a)=2 ≠ nrows(b)=1
assert_eq!(a.matmul(&b), Err(MatrixError::ShapeMismatch));
```

**参见**：[solve](#solve)

---

### `determinant`

**签名**：`pub fn determinant(&self) -> Result<D::Element, MatrixError>`

**功能**：使用 **Bareiss 无分数算法**（带部分选主元）计算方阵行列式。

**算法要点**：
- 通过消元将矩阵化为上三角，每步除以前一步的主元（Bareiss 保证整除）
- 部分选主元：遇到零主元时向下搜索非零行并交换
- 若消元过程中整列全零，行列式为零
- 符号由行交换次数的奇偶性决定

**返回值**：
- `Ok(det)` — 行列式值
- `Err(MatrixError::ShapeMismatch)` — 矩阵非方阵

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// 2×2: det = 1*4 - 2*3 = -2
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(3), Integer::from(4)]],
    d,
);
assert_eq!(a.determinant().unwrap(), Integer::from(-2));
```

**奇异矩阵**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// 第二行是第一行的 2 倍
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(2), Integer::from(4)]],
    d,
);
assert_eq!(a.determinant().unwrap(), Integer::from(0));
```

**参见**：[inverse](#inverse)、[rank](#rank)

---

### `rank`

**签名**：`pub fn rank(&self) -> usize`

**功能**：通过无分数高斯消元计算矩阵的秩。克隆矩阵后执行 `row_echelon`，返回非零行数。

**返回值**：`usize`，矩阵的秩。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// 满秩 2×2
let id = Matrix::identity(2, d);
assert_eq!(id.rank(), 2);

// 秩亏缺：第二行 = 2 × 第一行
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2), Integer::from(3)],
         vec![Integer::from(2), Integer::from(4), Integer::from(6)],
         vec![Integer::from(1), Integer::from(1), Integer::from(1)]],
    d,
);
assert_eq!(a.rank(), 2);
```

**参见**：[determinant](#determinant)、[row_echelon](#row_echelon)

---

### `inverse`

**签名**：`pub fn inverse(&self) -> Result<Matrix<D>, MatrixError>`

**功能**：计算方阵的逆矩阵 $A^{-1}$。通过逐列求解 $A x_j = e_j$（$e_j$ 为单位矩阵第 $j$ 列）实现。

**返回值**：
- `Ok(Matrix)` — 逆矩阵，满足 $A \cdot A^{-1} = I$
- `Err(MatrixError::ShapeMismatch)` — 矩阵非方阵
- `Err(MatrixError::Inconsistent)` — 矩阵奇异（无解）
- `Err(MatrixError::Underdetermined { .. })` — 奇异导致欠定
- `Err(MatrixError::ResultNotInDomain)` — 逆矩阵元素不在域中

**⚠️ 注意**：在整数域 `IntegerDomain` 上，仅单模矩阵（$\det = \pm 1$）有整数逆。一般矩阵在有理数域 `RationalDomain` 上求逆。

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// 单模矩阵：det = 1*1 - 0*2 = 1
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(0), Integer::from(1)]],
    d,
);
let inv = a.inverse().unwrap();
assert_eq!(inv[(0, 0)], Integer::from(1));
assert_eq!(inv[(0, 1)], Integer::from(-2));
assert_eq!(inv[(1, 0)], Integer::from(0));
assert_eq!(inv[(1, 1)], Integer::from(1));

// 验证 A * A^{-1} = I
let prod = a.matmul(&inv).unwrap();
assert_eq!(prod, Matrix::identity(2, IntegerDomain));
```

**奇异矩阵错误**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(2)],
         vec![Integer::from(2), Integer::from(4)]],
    d,
);
assert!(a.inverse().is_err());
```

**参见**：[solve](#solve)、[determinant](#determinant)

---

### `augment`

**签名**：`pub fn augment(&self, other: &Matrix<D>) -> Result<Matrix<D>, MatrixError>`

**功能**：水平拼接两个矩阵，得到增广矩阵 $[A \mid B]$。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `other` | `&Matrix<D>` | 要拼接在右侧的矩阵 |

**返回值**：
- `Ok(Matrix)` — 增广矩阵，形状为 `(nrows, ncols_a + ncols_b)`
- `Err(MatrixError::ShapeMismatch)` — 两矩阵行数不一致

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(3), Integer::from(2)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = Matrix::new(2, 1, vec![Integer::from(5), Integer::from(2)], d);
let aug = a.augment(&b).unwrap();
assert_eq!(aug.nrows(), 2);
assert_eq!(aug.ncols(), 3);
assert_eq!(aug[(0, 0)], Integer::from(3));
assert_eq!(aug[(0, 2)], Integer::from(5));
```

**参见**：[solve](#solve)

---

### `solve`

**签名**：`pub fn solve(&self, b: &[D::Element]) -> Result<Vec<D::Element>, MatrixError>`

**功能**：求解线性方程组 $Ax = b$。使用无分数高斯消元 + 回代。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `b` | `&[D::Element]` | 右端向量，长度必须等于 `self.nrows` |

**返回值**：
- `Ok(Vec<D::Element>)` — 解向量，长度等于 `self.ncols`
- `Err(MatrixError::ShapeMismatch)` — `b.len() != self.nrows`
- `Err(MatrixError::Inconsistent)` — 方程组无解
- `Err(MatrixError::Underdetermined { rank })` — 无穷多解，`rank < ncols`
- `Err(MatrixError::ResultNotInDomain)` — 解不在域中（如整数除法失败）

**算法步骤**：
1. 构造增广矩阵 $[A \mid b]$
2. 对前 `ncols` 列执行无分数行阶梯化（`row_echelon`）
3. 检查秩以下行的常数项是否为零（一致性检查）
4. 执行无分数回代（`back_substitution`）
5. 若秩 < 变量数，返回 `Underdetermined`
6. 最终除以主元得到解

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;

// 2x + y = 4
// x  + y = 3  → x = 1, y = 2
let a = Matrix::from_rows(
    vec![vec![Integer::from(2), Integer::from(1)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(4), Integer::from(3)];
let x = a.solve(&b).unwrap();
assert_eq!(x, vec![Integer::from(1), Integer::from(2)]);
```

**3×3 系统**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
// x + y + z = 6
// 2x - y + z = 3
// x + 2y - z = 2  → x=1, y=2, z=3
let a = Matrix::from_rows(
    vec![
        vec![Integer::from(1), Integer::from(1), Integer::from(1)],
        vec![Integer::from(2), Integer::from(-1), Integer::from(1)],
        vec![Integer::from(1), Integer::from(2), Integer::from(-1)],
    ],
    d,
);
let b = vec![Integer::from(6), Integer::from(3), Integer::from(2)];
let x = a.solve(&b).unwrap();
assert_eq!(x, vec![Integer::from(1), Integer::from(2), Integer::from(3)]);
```

**不相容系统**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(1)],
         vec![Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(1), Integer::from(2)];
assert_eq!(a.solve(&b), Err(MatrixError::Inconsistent));
```

**欠定系统**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::{Matrix, MatrixError};

let d = IntegerDomain;
let a = Matrix::from_rows(
    vec![vec![Integer::from(1), Integer::from(1), Integer::from(1)]],
    d,
);
let b = vec![Integer::from(3)];
assert!(matches!(a.solve(&b), Err(MatrixError::Underdetermined { .. })));
```

**参见**：[inverse](#inverse)、[augment](#augment)

---

## 内部消元方法

以下两个方法为 `pub`，通常不需要直接调用，但在某些高级场景下（如自定义消元管线）可能有用。

### `row_echelon`

**签名**：`pub fn row_echelon(&mut self, max_col: usize) -> usize`

**功能**：对前 `max_col` 列执行无分数高斯消元，将矩阵化为行阶梯形。就地修改矩阵。返回秩（非零行数）。

**参数**：

| 参数 | 类型 | 说明 |
|---|---|---|
| `max_col` | `usize` | 消元涉及的最大列数（不含），自动钳制到 `ncols` |

**返回值**：`usize` — 矩阵的秩。

**算法细节**：
- 每列寻找非零主元，必要时交换行
- 从主元行提取 GCD 内容以防止系数增长
- 消元采用 gcd 约化的无分数缩放：设 $g = \gcd(\text{pivot}, a_{kj})$，被消行乘以 $\text{pivot} / g$，主元行乘以 $a_{kj} / g$，相减后不含分数

**示例**：

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

let d = IntegerDomain;
let mut a = Matrix::from_rows(
    vec![
        vec![Integer::from(1), Integer::from(2), Integer::from(3)],
        vec![Integer::from(2), Integer::from(4), Integer::from(6)],
        vec![Integer::from(0), Integer::from(1), Integer::from(1)],
    ],
    d,
);
let rank = a.row_echelon(3);
assert_eq!(rank, 2);
```

**参见**：[back_substitution](#back_substitution)、[rank](#rank)

---

### `back_substitution`

**签名**：`pub fn back_substitution(&mut self, max_col: usize)`

**功能**：对已处于行阶梯形的矩阵执行无分数回代。就地修改矩阵的前 `max_col` 列。

**算法细节**：
- 从最后一行开始向上处理
- 对每个主元行，先提取 GCD 内容
- 对主元上方的行，用同样的 gcd 约化缩放（被消行乘以 $\text{pivot} / g$、主元行乘以 $a_{kj} / g$，其中 $g = \gcd(\text{pivot}, a_{kj})$）后消去主元列

**参见**：[row_echelon](#row_echelon)、[solve](#solve)

---

## 完整示例

以下综合示例演示矩阵的基本工作流：构造、运算、求解。

```rust
use ocas_domain::{Integer, IntegerDomain};
use ocas_poly::matrix::Matrix;

fn main() {
    let d = IntegerDomain;

    // 构造矩阵 A
    let a = Matrix::from_rows(
        vec![
            vec![Integer::from(2), Integer::from(1), Integer::from(-1)],
            vec![Integer::from(-3), Integer::from(-1), Integer::from(2)],
            vec![Integer::from(-2), Integer::from(1), Integer::from(2)],
        ],
        d,
    );

    // 行列式
    let det = a.determinant().unwrap();
    assert_eq!(det, Integer::from(-1)); // det = -1（单模）

    // 秩
    assert_eq!(a.rank(), 3);

    // 迹
    let tr = a.trace().unwrap();
    assert_eq!(tr, Integer::from(3)); // 2 + (-1) + 2

    // 转置
    let at = a.transpose();
    assert_eq!(at[(0, 1)], Integer::from(-3));

    // 求解 Ax = b
    let b = vec![Integer::from(8), Integer::from(-11), Integer::from(-3)];
    let x = a.solve(&b).unwrap();
    assert_eq!(x, vec![Integer::from(2), Integer::from(3), Integer::from(-1)]);

    // 逆矩阵（det = -1，整数逆存在）
    let inv = a.inverse().unwrap();
    // 验证 A * A^{-1} = I
    let prod = a.matmul(&inv).unwrap();
    assert_eq!(prod, Matrix::identity(3, IntegerDomain));

    // 矩阵乘法
    let c = a.matmul(&a).unwrap();
    assert_eq!(c.nrows(), 3);
    assert_eq!(c.ncols(), 3);
}
```

---

## 参见

- [表达式系统](./rust-expressions.md) — 表达式树构造与操作
- [系数域](./rust-domains.md) — `EuclideanDomain` trait 与域实现
- [多项式](./rust-polynomials.md) — 多项式类型，使用相同的域 trait
- [线性代数基础](../math/linear-algebra.md) — Bareiss 算法与高斯消元的数学原理
