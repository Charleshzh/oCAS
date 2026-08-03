# 快速上手

本章帮助你在 5 分钟内完成 oCAS 的安装与首次运行。oCAS 0.24 支持 **Rust**、**Python** 和 **C/C++** 三种语言接口。

> **系统要求**
>
> | 语言 | 最低版本 | 备注 |
> |---|---|---|
> | Rust | 1.97+ (edition 2024) | 推荐通过 [rustup](https://rustup.rs/) 安装 |
> | Python | 3.10+ | 使用 `pip install` 或从源码构建 |
> | C/C++ | C11 / C++17 | 需要 `cargo` 构建 C 库 |

---

## Rust

### 安装

在你的 `Cargo.toml` 中添加 oCAS：

```toml
[dependencies]
ocas = "0.24"
```

如需可选后端，在 `features` 中开启：

```toml
[dependencies]
ocas = { version = "0.24", features = ["gmp", "jit"] }
```

常用 feature 一览：

| Feature | 说明 |
|---|---|
| `gmp` | 使用 GMP 任意精度算术（替代纯 Rust `num-bigint`） |
| `mpfr` | 使用 MPFR 球算术（实数区间算术） |
| `jit` | Cranelift JIT 编译求值器 |
| `simd` | SIMD 向量化批量求值 |
| `egg` | e-graph 等式饱和引擎 |
| `ntt` | NTT 快速多项式乘法（𝔽_p） |
| `mimalloc` | 使用 mimalloc 分配器 |
| `system-libs` | 链接系统预装的 GMP/MPFR（Windows/MinGW 必需） |

完整 feature 列表见 [Rust API 参考](./api/rust.md)。

### 验证安装

创建一个最小项目来确认一切正常：

```bash
cargo new ocas-test && cd ocas-test
```

在 `Cargo.toml` 中添加依赖：

```toml
[dependencies]
ocas = "0.24"
```

在 `src/main.rs` 中写入：

```rust
use ocas::prelude::*;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, "x^2 + 2*x + 1").unwrap();
    let d = diff(&ctx, expr, Symbol::new("x"));
    println!("d/dx(x² + 2x + 1) = {}", d);
    // 输出：d/dx(x² + 2x + 1) = 2 + 2*x
}
```

运行：

```bash
cargo run
```

若输出 `2 + 2*x`，说明安装成功。

### 基本用法

#### 符号微分与化简

```rust
use ocas::prelude::*;
use ocas::ocas_rewrite::rules::default_rules;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);

    // 解析表达式
    let expr = parse(&ctx, "x^2").unwrap();

    // 微分
    let d = diff(&ctx, expr, Symbol::new("x"));
    println!("微分: {}", d);
    // 输出：2*x

    // 化简（默认规则集：恒等元移除、常量折叠）
    let e2 = parse(&ctx, "x + 0 + y*0 + z*1").unwrap();
    let rules = default_rules(&ctx, &());
    let s = simplify(&ctx, e2, &rules, 100);
    println!("化简: {}", s);
    // 输出：x + z
}
```

#### 多项式运算

```rust
use ocas::prelude::*;
use ocas_poly::dense::DenseUnivariatePolynomial;
use ocas_domain::rational::RationalDomain;
use ocas_domain::Rational;

fn main() {
    let domain = RationalDomain;
    // 1 + 2x + x²
    let p = DenseUnivariatePolynomial::from_coeffs(
        domain,
        vec![Rational::new(1, 1), Rational::new(2, 1), Rational::new(1, 1)],
    );
    // 1 + x
    let q = DenseUnivariatePolynomial::from_coeffs(
        domain,
        vec![Rational::new(1, 1), Rational::new(1, 1)],
    );

    let (quot, rem) = p.div_rem(&q).unwrap();
    // DenseUnivariatePolynomial 未实现 Display，直接断言系数
    assert_eq!(quot.coeffs(), &[Rational::new(1, 1), Rational::new(1, 1)]); // 1 + x
    assert_eq!(rem.coeffs(), &[]); // 余 0
}
```

#### 求解线性方程组

```rust
use ocas_calc::solve::solve_linear_rational;

fn main() {
    // 求解 2x + y = 5, x - y = 1
    let a = vec![vec![2, 1], vec![1, -1]];
    let b = vec![5, 1];

    let solution = solve_linear_rational(&a, &b).unwrap();
    println!("x = {}/{}, y = {}/{}",
        solution[0].0, solution[0].1,
        solution[1].0, solution[1].1);
    // 输出：x = 2/1, y = 1/1
}
```

#### 数论

```rust
use ocas_domain::number_theory::{
    is_prime, factor::factor_integer, functions::euler_phi,
};

fn main() {
    println!("is_prime(97) = {}", is_prime(&97.into()));
    // 输出：is_prime(97) = true

    let factors = factor_integer(&360.into());
    println!("360 = {}", factors.iter()
        .map(|(p, e)| format!("{}^{}", p, e))
        .collect::<Vec<_>>()
        .join(" × "));
    // 输出：360 = 2^3 × 3^2 × 5^1

    println!("φ(360) = {}", euler_phi(&360.into()));
    // 输出：φ(360) = 96
}
```

更多 Rust API 详见：
- [表达式系统](./api/rust-expressions.md) — Arena、Symbol、Atom、模式匹配
- [系数域](./api/rust-domains.md) — Integer、Rational、FiniteField 等
- [多项式](./api/rust-polynomials.md) — 一元/多元多项式、单项式序
- [矩阵](./api/rust-matrix.md) — 行列式、求解、乘法
- [微积分](./api/rust-calculus.md) — diff、integrate、taylor
- [求解器](./api/rust-solvers.md) — 线性方程、丢番图、ODE
- [数论](./api/rust-ntheory.md) — 素性、分解、离散对数

---

## Python

### 安装

```bash
pip install ocas
```

> 要求 Python ≥ 3.10（使用 PyO3 abi3 稳定 ABI）。

### 验证安装

```bash
python -c "import ocas; e = ocas.Expression('x^2'); print(e.diff('x'))"
```

若输出 `2*x`，说明安装成功。如需在 Python 交互环境中验证：

```python
>>> import ocas
>>> e = ocas.Expression("x^2 + 2*x + 1")
>>> e.diff("x")
2 + 2*x
>>> e.simplify()
1 + 2*x + x^2
```

### 基本用法

#### 符号微分与化简

```python
import ocas

e = ocas.Expression("x^2")

# 微分
d = e.diff("x")
print(f"微分: {d}")
# 输出：微分: 2*x

# 化简（默认规则集：恒等元移除、常量折叠）
e2 = ocas.Expression("x + 0 + y*0 + z*1")
s = e2.simplify()
print(f"化简: {s}")
# 输出：化简: x + z
```

#### 多项式运算

```python
import ocas

# 1 + 2x + x²
p = ocas.Polynomial([1, 2, 1])
print(f"次数: {p.degree()}")         # 次数: 2
print(f"p(2) = {p.eval(2)}")         # p(2) = 9

# 多项式 GCD
q = ocas.Polynomial([1, 1])          # 1 + x
g = p.gcd(q)
print(f"gcd = {g}")                   # gcd = 1 + x

# 因式分解
factors = p.factor()
print(f"因式分解: {factors}")         # 因式分解: [(1 + x, 2)]
```

#### 有限域运算

```python
import ocas

gf5 = ocas.FiniteField(5)
q = ocas.Polynomial([1, 2, 1], domain=gf5)
print(q.eval(3))  # (1 + 6 + 9) mod 5 = 16 mod 5 = 1
```

#### 矩阵运算

```python
import ocas

m = ocas.Matrix([[1, 2], [3, 4]])
print(f"行列式: {m.determinant()}")       # 行列式: -2
print(f"转置: {(m.transpose()).rows()}")   # 转置: [[1, 3], [2, 4]]
print(f"m × m: {(m @ m).rows()}")          # m × m: [[7, 10], [15, 22]]
```

#### 数论

```python
import ocas

print(ocas.isprime(97))               # True
print(ocas.factorint(360))             # [("2", 3), ("3", 2), ("5", 1)]
print(ocas.totient(360))               # 96
print(ocas.nextprime(100))             # 101
print(ocas.mobius(30))                 # -1
print(ocas.divisor_count(360))         # 24
```

#### ODE 求解

```python
import ocas

e = ocas.Expression("Derivative(y(x), x) + y(x)")  # y' + y = 0

# 分类可用的求解方法
methods = ocas.classify_ode(e, "y", "x")
print(methods)

# 求解 y' + y = 0
sol = ocas.dsolve(e, "y", "x")
print(sol)  # y = C1*exp(-x)
```

更多 Python API 详见 [Python API 参考](./api/python.md)。

---

## C/C++

### 构建 C 库

oCAS 的 C 库通过 `cargo` 构建。从源码仓库克隆并编译：

```bash
git clone https://github.com/charleshzh/ocas.git
cd ocas
cargo build -p ocas-c --release
```

构建产物位于 `target/release/`：

| 文件 | 说明 |
|---|---|
| `libocas_c.so`（Linux）/ `ocas_c.dll`（Windows）/ `libocas_c.dylib`（macOS） | 动态库 |
| `libocas_c.a`（Linux/macOS）/ `ocas_c.lib`（Windows） | 静态库 |

头文件 `ocas.h` 位于 `ocas-c/include/` 目录。C++ 用户还可以使用同目录下的 `ocas.hpp`（RAII 包装，覆盖表达式、数论和张量域）。

> **Windows 用户**：推荐使用 MSYS2 MINGW64 环境。详见 [在 Windows 上构建](./build-windows.md)。

### 验证安装

创建一个最小测试文件 `test_ocas.c`：

```c
#include <stdio.h>
#include "ocas.h"

int main(void) {
    struct ocas_OcasExpr *e = ocas_expr_parse("x^2 + 2*x + 1", NULL);
    if (!e) {
        fprintf(stderr, "解析失败\n");
        return 1;
    }

    struct ocas_OcasExpr *d = ocas_expr_diff(e, "x", NULL);
    char *s = ocas_expr_to_string(d, NULL);
    printf("d/dx(x² + 2x + 1) = %s\n", s);
    /* 输出：d/dx(x² + 2x + 1) = 2 + 2*x */

    ocas_string_free(s);
    ocas_expr_free(d);
    ocas_expr_free(e);
    return 0;
}
```

编译并运行（以 Linux 为例，调整路径以匹配你的环境）：

```bash
# 假设 oCAS 已构建在 /path/to/ocas
gcc test_ocas.c -I/path/to/ocas/ocas-c/include \
    -L/path/to/ocas/target/release -locas_c -lm -ldl -lpthread -o test_ocas
LD_LIBRARY_PATH=/path/to/ocas/target/release ./test_ocas
```

若输出 `2 + 2*x`，说明安装成功。

### 基本用法

#### 符号微分

```c
#include <stdio.h>
#include "ocas.h"

int main(void) {
    struct ocas_OcasExpr *e = ocas_expr_parse("x^2", NULL);
    struct ocas_OcasExpr *d = ocas_expr_diff(e, "x", NULL);

    char *s = ocas_expr_to_string(d, NULL);
    printf("微分: %s\n", s);
    /* 输出：微分: 2*x */

    ocas_string_free(s);
    ocas_expr_free(d);
    ocas_expr_free(e);
    return 0;
}
```

#### 数论

```c
#include <stdio.h>
#include "ocas.h"

int main(void) {
    char *result;

    /* 素性测试 */
    int prime = ocas_ntheory_isprime("97", NULL);
    printf("is_prime(97) = %s\n", prime ? "true" : "false");
    /* 输出：is_prime(97) = true */

    /* 整数分解 */
    result = ocas_ntheory_factorint("360", NULL);
    printf("360 = %s\n", result);
    /* 输出：360 = 2:3,3:2,5:1 */
    ocas_string_free(result);

    /* Euler φ 函数 */
    result = ocas_ntheory_totient("360", NULL);
    printf("φ(360) = %s\n", result);
    /* 输出：φ(360) = 96 */
    ocas_string_free(result);

    return 0;
}
```

#### C++ RAII 包装

如果你使用 C++，可以利用 `ocas.hpp` 中的 RAII 包装避免手动管理内存：

```cpp
#include <iostream>
#include "ocas.h"
#include "ocas.hpp"

int main() {
    ocas::Expression e("x^2 + 2*x + 1");
    auto d = e.diff("x");
    std::cout << "d/dx = " << d.to_string() << std::endl;
    // 输出：d/dx = 2 + 2*x

    // 数论（无需手动 free）
    std::cout << "is_prime(97) = " << ocas::ntheory::isprime("97") << std::endl;
    std::cout << "360 = " << ocas::ntheory::factorint("360") << std::endl;
    std::cout << "φ(360) = " << ocas::ntheory::totient("360") << std::endl;

    return 0;
}
```

更多 C/C++ API 详见 [C/C++ API 参考](./api/c.md)。

---

## 常见问题

### Rust 相关

**Q：编译很慢怎么办？**

oCAS 包含大量泛型代码，首次编译可能较慢。建议：
- 开发时使用 `cargo build`（dev profile，已对依赖开 `opt-level = 3`）
- 发布时使用 `cargo build --release`（启用 LTO）
- 利用 `sccache` 或 `cargo-nextest` 加速增量编译

**Q：`gmp` feature 报错找不到 GMP 库？**

在 Linux/macOS 上，确保已安装 GMP 开发包：
```bash
# Ubuntu/Debian
sudo apt install libgmp-dev

# macOS (Homebrew)
brew install gmp
```

在 Windows（MSYS2 MINGW64）上：
```bash
pacman -S mingw-w64-x86_64-gmp
```

或使用 `system-libs` feature 链接系统预装库。详见 [在 Windows 上构建](./build-windows.md)。

**Q：`flint` feature 在 Windows 上不可用？**

上游 `flint3-sys` 仅支持 POSIX 平台，Windows 上暂不支持 `flint` feature。

### Python 相关

**Q：`pip install ocas` 失败？**

确认 Python 版本 ≥ 3.10：
```bash
python --version
```

如需从源码构建：
```bash
pip install maturin
maturin develop --release -m ocas-py/Cargo.toml
```

**Q：Python 绑定中的函数名为什么有 `py_` 前缀？**

部分 Gröbner 函数（如 `ocas.py_groebner_basis`）保留了 `py_` 前缀，这是因为 PyO3 绑定尚未添加 `#[pyo3(name=...)]` 注解。实际调用时使用带前缀的名称即可。

### C/C++ 相关

**Q：链接时报 undefined reference？**

确保链接顺序正确——`libocas_c` 放在源文件之后：
```bash
gcc test.c -locas_c -lm -ldl -lpthread
```

并确保 `LD_LIBRARY_PATH`（Linux）或 `PATH`（Windows）包含库文件所在目录。

**Q：`ocas.h` 中缺少某些函数？**

仓库中附带的 `ocas.h` 可能不是最新的。构建时 cbindgen 会在 `target/release/build/ocas-c-*/out/` 下生成最新的头文件。以该文件为准。

**Q：C++ 包装 `ocas.hpp` 覆盖了哪些功能？**

`ocas.hpp` 目前覆盖三个域：**表达式**（`ocas::Expression`）、**数论**（`ocas::ntheory::*`）和**张量**（`ocas::tensor::*`）。其余域（多项式、Gröbner、ODE 等）需直接使用 C API。

### 通用问题

**Q：数值结果不精确？**

默认使用纯 Rust `num-bigint` 有理数算术，完全精确。若需高精度浮点（区间算术），启用 `mpfr` feature：
```toml
ocas = { version = "0.24", features = ["mpfr"] }
```

**Q：如何查看 oCAS 版本？**

```bash
cargo metadata --format-version=1 | jq '.packages[] | select(.name=="ocas") | .version'
```

---

## 下一步

- [Rust API 参考](./api/rust.md) — 完整的 Rust 模块、类型与函数文档
- [Python API 参考](./api/python.md) — 所有 Python 类与函数的详细说明
- [C/C++ API 参考](./api/c.md) — C 导出函数、C++ RAII 包装、内存管理
- [架构](./architecture.md) — oCAS 的整体设计与模块关系
- [Gröbner 基实现](./algorithms/groebner.md) — Buchberger / F4 / F5 算法详解
- [ODE 求解](./algorithms/ode-solving.md) — 常微分方程分类与求解方法
