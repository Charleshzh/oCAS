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
