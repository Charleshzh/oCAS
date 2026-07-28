# C/C++ API

`ocas-c` crate 提供稳定的 C ABI（由 `cbindgen` 生成），覆盖表达式生命周期、微积分与化简，并在 `ocas-c/include/ocas.hpp` 提供 C++ RAII 包装。

## 构建

```bash
cargo build -p ocas-c --release
```

共享库与 `ocas.h` / `ocas.hpp` 头文件位于 `ocas-c/include/`。

## C 示例

```c
#include <ocas.h>

int main(void) {
    ocas_error err;
    ocas_expr* e = ocas_expr_parse("x^2 + 2*x + 1", &err);
    ocas_expr* d = ocas_expr_diff(e, "x", &err);

    char* s = ocas_expr_to_string(d, &err);
    printf("derivative: %s\n", s);   /* 2*x + 2 */

    ocas_string_free(s);
    ocas_expr_free(d);
    ocas_expr_free(e);
    return 0;
}
```

## C++ RAII

```cpp
#include <ocas.hpp>

int main() {
    ocas::Expression e("x^2 + 2*x + 1");
    auto d = e.diff("x");
    std::cout << d.to_string() << std::endl;   // 2*x + 2
    return 0;   // 自动清理
}
```

C++ 包装将 oCAS 错误转换为 `ocas::Error` 异常，并通过 RAII 管理 arena 后端表达式，无需手动调用 `free`。

## 多项式 API

自 0.11.1 起，`ocas-c` 将二元多项式对象暴露为不透明句柄，支持 $\mathbb{Z}$ 和 $\mathbb{F}_p$ 上的因式分解。

### 整数多项式（`OcasPolyZ`）

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err;
    // 从字符串创建二元整数多项式。
    OcasPolyZ* p = ocas_poly_z_create("x^2 + y + 1", &err);

    // 查询总次数。
    printf("degree: %zu\n", ocas_poly_z_degree(p));

    // 因式分解为不可约因子。
    OcasPolyFactorArray factors;
    ocas_poly_z_factor(p, &factors, &err);

    for (size_t i = 0; i < factors.len; i++) {
        OcasPolyZ* fi = (OcasPolyZ*)factors.factors[i].poly;
        char* s = ocas_poly_z_to_string(fi, &err);
        printf("  factor %zu: %s (mult %zu)\n", i, s,
               factors.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_z_free(fi);
    }
    ocas_poly_factor_array_free(&factors);
    ocas_poly_z_free(p);
    return 0;
}
```

### 有限域多项式（`OcasPolyFp`）

```c
// 在 F_5 上创建多项式并因式分解。
OcasPolyFp* p = ocas_poly_fp_create("x^2 + y + 1", "5", &err);

OcasPolyFactorArray factors;
ocas_poly_fp_factor(p, &factors, &err);

for (size_t i = 0; i < factors.len; i++) {
    OcasPolyFp* fi = (OcasPolyFp*)factors.factors[i].poly;
    char* s = ocas_poly_fp_to_string(fi, &err);
    printf("  factor %zu: %s\n", i, s);
    ocas_string_free(s);
    ocas_poly_fp_free(fi);
}
ocas_poly_factor_array_free(&factors);
ocas_poly_fp_free(p);
```

### 生命周期

| 函数 | 用途 |
|---|---|
| `ocas_poly_z_create` / `ocas_poly_fp_create` | 从字符串创建 |
| `ocas_poly_z_clone` / `ocas_poly_fp_clone` | 深拷贝 |
| `ocas_poly_z_degree` / `ocas_poly_fp_degree` | 总次数 |
| `ocas_poly_z_to_string` / `ocas_poly_fp_to_string` | 堆分配字符串（调用者释放） |
| `ocas_poly_z_factor` / `ocas_poly_fp_factor` | 因式分解为不可约因子 |
| `ocas_poly_z_free` / `ocas_poly_fp_free` | 释放句柄 |
| `ocas_poly_factor_array_free` | 释放因子数组 |

所有多项式函数均可安全调用（无需 `unsafe`）。传入 `NULL` 会设置错误码并返回 `NULL` / 错误。

## ODE API

自 0.20.1 起，`ocas-c` 以字符串进出方式暴露 ODE 求解。所有返回字符串
均为堆分配，须用 `ocas_string_free` 释放。

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // 只分类不求解。
    char *types = ocas_ode_classify(
        "Derivative(y(x), x) - y(x)", "y", "x", &err);
    printf("methods: %s\n", types);   /* LinearFirst,Separable,PowerSeries */

    // 求解（hint = NULL 自动分类）。
    char *sol = ocas_ode_dsolve(
        "Derivative(y(x), x) - y(x)", "y", "x", NULL, &err);
    printf("solution: %s\n", sol);    /* y = C1*exp(x) */

    // Laplace 初值问题：y(0) = 1。
    char *ivp = ocas_ode_dsolve_ivp(
        "Derivative(y(x), x) - y(x)", "y", "x", "1", NULL, &err);
    printf("ivp: %s\n", ivp);          /* y = exp(x) */

    ocas_string_free(types);
    ocas_string_free(sol);
    ocas_string_free(ivp);
    return 0;
}
```

| 函数 | 用途 |
|---|---|
| `ocas_ode_classify` | 逗号分隔的可用方法名 |
| `ocas_ode_dsolve` | 符号解字符串（`"y = ..."` 或 `"unsolved"`） |
| `ocas_ode_dsolve_ivp` | Laplace 变换求显式 IVP 解 |

## 数值积分 API（Vegas）

自 0.18.0 起，`ocas-c` 通过不透明句柄提供蒙特卡洛积分。

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // 创建 2 维 Vegas 积分器（n_samples=10000, iterations=10）。
    ocas_OcasVegas* v = ocas_vegas_create(2, 10000, 10, &err);

    // 使用回调积分。
    ocas_OcasIntegrateResult result;
    ocas_integrate_1d(my_fn, NULL, 0.0, 1.0, 10000, 10, &result);
    printf("integral = %f ± %f\n", result.integral, result.error);

    ocas_vegas_free(v);
    return 0;
}
```

| 函数 | 用途 |
|---|---|
| `ocas_vegas_create(n_dims, n_samples, iterations, &err)` | 创建积分器 |
| `ocas_vegas_integrate(vegas, fn, user_data, &err)` | 运行积分 |
| `ocas_vegas_result(vegas)` | 获取 `OcasIntegrateResult` |
| `ocas_integrate_1d(fn, user_data, a, b, n_samples, iterations, &result)` | 一维便捷函数 |

## 双数 API（HyperDual）

自 0.18.1 起，`ocas-c` 通过不透明句柄提供前向自动微分。系数为字符串
（`"num"` 或 `"num/den"`）。

```c
#include <ocas.h>

int main(void) {
    int err = 0;
    ocas_OcasDualShape* shape = ocas_dual_shape_new(2, &err);

    // x = 1 + ε₁, y = 2 + ε₂
    ocas_OcasHyperDual* x = ocas_dual_variable(shape, 0, "1", &err);
    ocas_OcasHyperDual* y = ocas_dual_variable(shape, 1, "2", &err);

    // f = x * y
    ocas_OcasHyperDual* f = ocas_dual_mul(x, y, &err);
    char* val = ocas_dual_value(f, &err);    /* "2" */
    char* dx  = ocas_dual_deriv(f, 0, &err); /* "2" */

    ocas_string_free(val);
    ocas_string_free(dx);
    ocas_hyperdual_free(f);
    ocas_hyperdual_free(y);
    ocas_hyperdual_free(x);
    ocas_dual_shape_free(shape);
    return 0;
}
```

| 函数 | 用途 |
|---|---|
| `ocas_dual_shape_new(n_vars, &err)` | 一阶形状 |
| `ocas_dual_variable(shape, i, coeff, &err)` | 带单位 ε 的变量 |
| `ocas_dual_constant(shape, coeff, &err)` | 常量（无 ε） |
| `ocas_dual_value(hd, &err)` | 标量值字符串 |
| `ocas_dual_deriv(hd, i, &err)` | 第 i 阶导数字符串 |
| `ocas_dual_add/sub/mul/div(a, b, &err)` | 算术运算 |

## 张量 API

自 0.18.1 起，`ocas-c` 提供张量创建、缩并与对称化。指标标签和位置以
数组传入。

```c
#include <ocas.h>

int main(void) {
    int err = 0;

    // 创建 A^μ（1 阶，上指标 "mu"）。
    const char* labels[] = {"mu"};
    int positions[] = {1};  /* 1 = upper */
    ocas_OcasTensor* A = ocas_tensor_create("A", labels, positions, 1, &err);

    int rank = ocas_tensor_rank(A, &err);  /* 1 */

    // 缩并两个张量（返回缩并结果）。
    ocas_OcasTensorContraction* c = ocas_tensor_contract(A, B, &err);
    /* 通过 c->scalar / c->product 访问标量或自由因子 */

    ocas_tensor_contraction_free(c);
    ocas_tensor_free(A);
    return 0;
}
```

| 函数 | 用途 |
|---|---|
| `ocas_tensor_create(name, labels, positions, rank, &err)` | 创建张量 |
| `ocas_tensor_rank(tensor, &err)` | 查询阶数 |
| `ocas_tensor_symmetry(tensor, &err)` | 查询对称性（0=None, 1=Sym, 2=Anti） |
| `ocas_tensor_symmetrise_sign(tensor, &err)` | 反对称化符号 |
| `ocas_tensor_contract(a, b, &err)` | 缩并匹配指标 |
| `ocas_tensor_to_string(tensor, &err)` | 字符串表示 |

## 代数数 API

自 0.17.1 起，`ocas-c` 通过不透明句柄提供代数数域运算。

```c
#include <ocas.h>

int main(void) {
    int err = 0;

    // Q(√2)：极小多项式 x^2 - 2
    ocas_OcasAlgebraicField* field =
        ocas_algebraic_field_create("x^2 - 2", &err);
    int deg = ocas_algebraic_field_degree(field, &err);  /* 2 */

    // 在 Q(√2) 上创建多项式并因式分解。
    ocas_OcasAlgebraicPoly* p =
        ocas_algebraic_poly_create(field, coeffs, n_coeffs, &err);
    ocas_OcasAlgebraicFactorArray factors;
    ocas_algebraic_poly_factor(p, &factors, &err);

    /* ... 遍历因子 ... */

    ocas_algebraic_factor_array_free(&factors);
    ocas_algebraic_poly_free(p);
    ocas_algebraic_field_free(field);
    return 0;
}
```

| 函数 | 用途 |
|---|---|
| `ocas_algebraic_field_create(min_poly, &err)` | 从极小多项式创建 |
| `ocas_algebraic_field_degree(field, &err)` | 扩张次数 |
| `ocas_algebraic_poly_create(field, coeffs, n, &err)` | 域上多项式 |
| `ocas_algebraic_poly_factor(poly, &factors, &err)` | 因式分解（Trager 算法） |
