# C/C++ API 参考

> **⚠️ 重要提示：头文件 `ocas.h` 已过时**
>
> 由 cbindgen 自动生成的 `ocas-c/include/ocas.h` 缺少以下 9 个函数：
> - `ocas_expr_integrate_heuristic`（表达式域）
> - 8 个 Gröbner 域导出（`ocas_groebner_basis`、`ocas_groebner_basis_free`、`ocas_groebner_basis_len`、`ocas_is_zero_dimensional`、`ocas_solve_polynomial_system`、`ocas_system_solution_count`、`ocas_system_solution_value`、`ocas_system_solution_free`）
>
> 本文档基于 `ocas-c/src/*.rs` 中的 **实际导出**，而非头文件。使用缺失函数时需自行声明函数原型。

> **⚠️ Gröbner 模块注意事项**
>
> Gröbner 模块使用**非标准** `-1` 错误码（而非 `OCAS_*` 常量），且 `OcasGroebnerBasis` 和 `OcasSystemSolution` 句柄**未标记** `#[repr(C)]`——它们在当前 Rust 版本中工作正常，但 ABI 稳定性不由编译器保证。

> **⚠️ C++ 包装覆盖范围**
>
> `ocas.hpp` 仅覆盖三个域（共 22 个包装函数）：
> - **表达式**：`ocas::Expression` 类（8 个成员：2 个构造函数 + 6 个方法）
> - **数论**：`ocas::ntheory::*`（11 个内联函数）
> - **张量**：`ocas::tensor::*`（3 个内联函数）
>
> 其余 69 个导出无 C++ RAII 包装，需直接调用 C API。

---

## 通用约定

### 错误处理

所有 fallible 函数遵循以下模式：

```c
int err = 0;
OcasExpr *e = ocas_expr_parse("x^2 + 1", &err);
if (e == NULL) {
    fprintf(stderr, "错误 %d: %s\n", err, ocas_error_last_message());
    ocas_error_clear();
    return 1;
}
```

- `err_out` 参数可为 `NULL`（忽略错误码）。
- 返回指针的函数：失败返回 `NULL`。
- 返回 `int` 的函数：失败返回非零值（具体含义见各函数文档）。
- 返回 `size_t` 的函数：失败返回 `0`。
- 返回 `int` 的数论函数使用特殊哨兵值 `-1`（`isprime`）或 `-2`（`jacobi`/`mobius`/`liouville`）表示错误。

### 内存管理

| 返回类型 | 释放函数 |
|---|---|
| `char*`（字符串） | `ocas_string_free()` |
| `OcasExpr*` | `ocas_expr_free()` |
| `OcasPolyZ*` | `ocas_poly_z_free()` |
| `OcasPolyFp*` | `ocas_poly_fp_free()` |
| `OcasPolyFactorArray` | `ocas_poly_factor_array_free()`（仅释放数组结构，不释放内部多项式句柄） |
| `OcasAlgebraicField*` | `ocas_algebraic_field_free()` |
| `OcasAlgebraicPoly*` | `ocas_algebraic_poly_free()` |
| `OcasAlgebraicFactorArray` | `ocas_algebraic_factor_array_free()`（同上，仅释放数组） |
| `OcasVegas*` | `ocas_vegas_free()` |
| `OcasTensor*` | `ocas_tensor_free()` |
| `OcasTensorContraction` | `ocas_tensor_contraction_free()`（释放数组/字符串，不释放内部张量句柄） |
| `OcasDualShape*` | `ocas_dual_shape_free()` |
| `OcasHyperDual*` | `ocas_hyperdual_free()` |
| `OcasGroebnerBasis*` | `ocas_groebner_basis_free()` |
| `OcasSystemSolution*` | `ocas_system_solution_free()` |

所有 `*_free` 函数对 `NULL` 参数安全（no-op）。

---

## 1. 错误/工具（3 个函数）

### `ocas_version`

**签名**：
```c
const char *ocas_version(void);
```

**功能**：返回 oCAS 版本字符串。

**返回值**：程序生命周期内有效的静态字符串指针。**调用方不得释放或修改**。

**示例**：
```c
printf("oCAS 版本: %s\n", ocas_version());
// 输出：oCAS 版本: 0.24.x
```

---

### `ocas_error_last_message`

**签名**：
```c
const char *ocas_error_last_message(void);
```

**功能**：返回调用线程上最后一次错误的消息，无错误时返回 `NULL`。

**返回值**：库拥有的字符串指针。**不得释放或修改**。在同一线程上调用任何可能设置错误的函数或 `ocas_error_clear()` 后失效。

**示例**：
```c
ocas_expr_parse("invalid !!!", NULL);
const char *msg = ocas_error_last_message();
if (msg) printf("错误: %s\n", msg);
ocas_error_clear();
```

---

### `ocas_error_clear`

**签名**：
```c
void ocas_error_clear(void);
```

**功能**：清除调用线程上的最后错误状态。

**示例**：
```c
ocas_error_clear(); // 重置错误状态
```

---

## 2. 表达式（12 个函数）

表达式是 oCAS 的核心类型。每个 `OcasExpr*` 句柄拥有私有 arena、一个 `AtomArena` 和一个规范化的根 `Atom`。句柄可移动；堆地址在句柄生命周期内不变。

### `ocas_expr_parse`

**签名**：
```c
OcasExpr *ocas_expr_parse(const char *input, int *err_out);
```

**功能**：将输入字符串解析为新的表达式句柄。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `input` | `const char*` | 以 null 结尾的表达式字符串，如 `"x^2 + 2*x + 1"` |
| `err_out` | `int*` | 可为 `NULL`；非空时写入错误码 |

**返回值**：成功返回 `OcasExpr*`；失败返回 `NULL`（`OCAS_ERROR_PARSE`）。

**内存**：调用方负责用 `ocas_expr_free()` 释放。

**示例**：
```c
int err = 0;
OcasExpr *e = ocas_expr_parse("sin(x)^2 + cos(x)^2", &err);
assert(err == OCAS_OK && e != NULL);
ocas_expr_free(e);
```

---

### `ocas_expr_free`

**签名**：
```c
void ocas_expr_free(OcasExpr *handle);
```

**功能**：释放表达式句柄。`NULL` 安全（no-op）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `OcasExpr*` | 由 `ocas_expr_parse` 或其他表达式函数返回的句柄 |

**示例**：
```c
ocas_expr_free(e);     // 释放
ocas_expr_free(NULL);  // no-op
```

---

### `ocas_expr_clone`

**签名**：
```c
OcasExpr *ocas_expr_clone(const OcasExpr *handle, int *err_out);
```

**功能**：克隆表达式到新 arena。返回独立副本。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `const OcasExpr*` | 非空表达式句柄 |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：新句柄或 `NULL`（`OCAS_ERROR_RUNTIME`）。

**示例**：
```c
OcasExpr *copy = ocas_expr_clone(original, NULL);
// copy 和 original 独立，互不影响
ocas_expr_free(copy);
```

---

### `ocas_expr_to_string`

**签名**：
```c
char *ocas_expr_to_string(const OcasExpr *handle, int *err_out);
```

**功能**：将表达式渲染为 null 结尾的 C 字符串。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `const OcasExpr*` | 非空表达式句柄 |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：堆分配的字符串（调用方用 `ocas_string_free()` 释放），或 `NULL`。

**示例**：
```c
char *s = ocas_expr_to_string(e, NULL);
printf("表达式: %s\n", s);
ocas_string_free(s);
```

---

### `ocas_string_free`

**签名**：
```c
void ocas_string_free(char *s);
```

**功能**：释放由 `ocas_expr_to_string` 或其他返回 `char*` 的函数分配的字符串。`NULL` 安全。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `s` | `char*` | 堆分配的字符串指针 |

---

### `ocas_expr_normalize`

**签名**：
```c
int ocas_expr_normalize(OcasExpr *handle, int *err_out);
```

**功能**：就地重新规范化表达式句柄。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `OcasExpr*` | 非空表达式句柄（**修改**） |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：`OCAS_OK`（0）成功；非零失败。

---

### `ocas_expr_diff`

**签名**：
```c
OcasExpr *ocas_expr_diff(const OcasExpr *handle, const char *var, int *err_out);
```

**功能**：对表达式关于变量 `var` 求导。返回新句柄。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `const OcasExpr*` | 非空表达式句柄 |
| `var` | `const char*` | 变量名，如 `"x"` |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：导数表达式句柄或 `NULL`。

**示例**：
```c
OcasExpr *e = ocas_expr_parse("x^3 + 2*x", NULL);
OcasExpr *de = ocas_expr_diff(e, "x", NULL);
char *s = ocas_expr_to_string(de, NULL);
printf("d/dx = %s\n", s);  // 输出：2 + (3*(x^2))
ocas_string_free(s);
ocas_expr_free(de);
ocas_expr_free(e);
```

---

### `ocas_expr_integrate`

**签名**：
```c
OcasExpr *ocas_expr_integrate(const OcasExpr *handle, const char *var, int *err_out);
```

**功能**：对表达式关于变量 `var` 积分。若无法解析求解，返回未求值形式 `Integral(expr, var)`。

**参数**：同 `ocas_expr_diff`。

**返回值**：积分表达式句柄或 `NULL`。

**示例**：
```c
OcasExpr *e = ocas_expr_parse("3*x^2", NULL);
OcasExpr *ie = ocas_expr_integrate(e, "x", NULL);
char *s = ocas_expr_to_string(ie, NULL);
printf("∫ = %s\n", s);  // 输出：3*(3^-1)*(x^3)
ocas_string_free(s);
ocas_expr_free(ie);
ocas_expr_free(e);
```

---

### `ocas_expr_integrate_heuristic`

**签名**：
```c
OcasExpr *ocas_expr_integrate_heuristic(const OcasExpr *handle, const char *var, int *err_out);
```

**功能**：使用启发式技术（分部积分、三角替换、Weierstrass 有理化、Euler 替换）积分。若无启发式成功，返回未求值形式。

> **⚠️ 注意**：此函数**未包含在 `ocas.h` 头文件中**。使用时需自行声明：
> ```c
> extern OcasExpr *ocas_expr_integrate_heuristic(const OcasExpr *handle,
>                                                 const char *var,
>                                                 int *err_out);
> ```

**参数**：同 `ocas_expr_diff`。

**返回值**：积分表达式句柄或 `NULL`。

---

### `ocas_expr_taylor`

**签名**：
```c
OcasExpr *ocas_expr_taylor(const OcasExpr *handle, const char *var,
                           const OcasExpr *point, uint32_t order, int *err_out);
```

**功能**：计算表达式在 `point` 处关于 `var` 的 `order` 阶 Taylor 展开。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `const OcasExpr*` | 非空表达式句柄 |
| `var` | `const char*` | 展开变量 |
| `point` | `const OcasExpr*` | 展开点（表达式句柄） |
| `order` | `uint32_t` | 展开阶数 |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：Taylor 展开句柄或 `NULL`。

**示例**：
```c
OcasExpr *f = ocas_expr_parse("exp(x)", NULL);
OcasExpr *zero = ocas_expr_parse("0", NULL);
OcasExpr *t = ocas_expr_taylor(f, "x", zero, 4, NULL);
char *s = ocas_expr_to_string(t, NULL);
printf("Taylor = %s\n", s);  // 输出：1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3)) + ((24^-1)*(x^4))
ocas_string_free(s);
ocas_expr_free(t);
ocas_expr_free(zero);
ocas_expr_free(f);
```

---

### `ocas_expr_simplify`

**签名**：
```c
OcasExpr *ocas_expr_simplify(const OcasExpr *handle, int *err_out);
```

**功能**：使用默认规则集化简表达式。返回新句柄。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `const OcasExpr*` | 非空表达式句柄 |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：化简后的表达式句柄或 `NULL`。

**示例**：
```c
OcasExpr *e = ocas_expr_parse("x + x", NULL);
OcasExpr *simplified = ocas_expr_simplify(e, NULL);
char *s = ocas_expr_to_string(simplified, NULL);
printf("化简 = %s\n", s);  // 输出：2*x
ocas_string_free(s);
ocas_expr_free(simplified);
ocas_expr_free(e);
```

---

### `ocas_expr_substitute`

**签名**：
```c
OcasExpr *ocas_expr_substitute(const OcasExpr *handle, const char *var,
                                const OcasExpr *replacement, int *err_out);
```

**功能**：将表达式中所有 `var` 的出现替换为 `replacement`。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `handle` | `const OcasExpr*` | 非空表达式句柄 |
| `var` | `const char*` | 被替换的变量名 |
| `replacement` | `const OcasExpr*` | 替换表达式（非空） |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：替换后的表达式句柄或 `NULL`。

**示例**：
```c
OcasExpr *e = ocas_expr_parse("x^2 + y", NULL);
OcasExpr *val = ocas_expr_parse("3", NULL);
OcasExpr *result = ocas_expr_substitute(e, "x", val, NULL);
char *s = ocas_expr_to_string(result, NULL);
printf("替换后 = %s\n", s);  // 输出：y + (3^2)
ocas_string_free(s);
ocas_expr_free(result);
ocas_expr_free(val);
ocas_expr_free(e);
```

---

### C++ 包装：`ocas::Expression`

`ocas.hpp` 提供 RAII 包装类，自动管理句柄生命周期：

```cpp
#include <ocas.h>
#include <ocas.hpp>

try {
    ocas::Expression e("x^2 + 1");
    ocas::Expression d = e.diff("x");
    std::cout << "d/dx = " << d.to_string() << std::endl;
    // 自动释放所有句柄
} catch (const ocas::Error& ex) {
    std::cerr << "错误: " << ex.what() << std::endl;
}
```

| C++ 方法 | 对应 C 函数 |
|---|---|
| `Expression(input)` | `ocas_expr_parse` |
| `Expression(other)` | `ocas_expr_clone` |
| `to_string()` | `ocas_expr_to_string` + `ocas_string_free` |
| `diff(var)` | `ocas_expr_diff` |
| `integrate(var)` | `ocas_expr_integrate` |
| `simplify()` | `ocas_expr_simplify` |
| `substitute(var, rep)` | `ocas_expr_substitute` |
| `raw()` | 获取原始 `OcasExpr*`（不转移所有权） |

---

## 3. 整数多项式（6 个函数）

处理双变量（`x`, `y`）整数多项式 `ℤ[x,y]`。从字符串解析，支持因式分解。

### `ocas_poly_z_create`

**签名**：
```c
OcasPolyZ *ocas_poly_z_create(const char *input, int *err);
```

**功能**：从字符串创建双变量整数多项式。输入可含变量 `x`、`y`，整数系数，加法、乘法和非负整数幂。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `input` | `const char*` | 多项式字符串，如 `"x^2*y + 3*x - 1"` |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasPolyZ*` 句柄或 `NULL`（`OCAS_ERROR_PARSE`）。

**内存**：调用方用 `ocas_poly_z_free()` 释放。

**示例**：
```c
int err = 0;
OcasPolyZ *p = ocas_poly_z_create("x^2 + y^2 - 1", &err);
assert(p != NULL);
ocas_poly_z_free(p);
```

---

### `ocas_poly_z_free`

**签名**：
```c
void ocas_poly_z_free(OcasPolyZ *poly);
```

**功能**：释放整数多项式句柄。`NULL` 安全。

---

### `ocas_poly_z_clone`

**签名**：
```c
OcasPolyZ *ocas_poly_z_clone(const OcasPolyZ *poly);
```

**功能**：克隆整数多项式。返回新句柄（调用方释放）。

**返回值**：新 `OcasPolyZ*` 或 `NULL`（输入为 `NULL` 或内存不足）。

---

### `ocas_poly_z_degree`

**签名**：
```c
size_t ocas_poly_z_degree(const OcasPolyZ *poly);
```

**功能**：返回多项式的**总次数**。零多项式返回 `0`。`NULL` 句柄返回 `0`。

---

### `ocas_poly_z_to_string`

**签名**：
```c
char *ocas_poly_z_to_string(const OcasPolyZ *poly, int *err);
```

**功能**：返回多项式的堆分配字符串表示。调用方用 `ocas_string_free()` 释放。

**返回值**：字符串或 `NULL`。

**示例**：
```c
char *s = ocas_poly_z_to_string(p, NULL);
printf("多项式: %s\n", s);
ocas_string_free(s);
```

---

### `ocas_poly_z_factor`

**签名**：
```c
int ocas_poly_z_factor(const OcasPolyZ *poly, OcasPolyFactorArray *out, int *err);
```

**功能**：对双变量整数多项式进行因式分解。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `poly` | `const OcasPolyZ*` | 非空多项式句柄 |
| `out` | `OcasPolyFactorArray*` | 输出结构体（成功时填充） |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OCAS_OK`（0）成功；非零失败（`out` 不变）。

**内存**：成功时，`out.factors[i].poly` 是 `void*`，需转换为 `OcasPolyZ*` 并用 `ocas_poly_z_free()` 释放。数组本身用 `ocas_poly_factor_array_free()` 释放。

**示例**：
```c
OcasPolyZ *p = ocas_poly_z_create("x^2 - y^2", NULL);
OcasPolyFactorArray arr = {0};
if (ocas_poly_z_factor(p, &arr, NULL) == OCAS_OK) {
    for (size_t i = 0; i < arr.len; i++) {
        OcasPolyZ *f = (OcasPolyZ *)arr.factors[i].poly;
        char *s = ocas_poly_z_to_string(f, NULL);
        printf("因子: %s (重数 %zu)\n", s, arr.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_z_free(f);
    }
    ocas_poly_factor_array_free(&arr);
}
ocas_poly_z_free(p);
```

---

## 4. 有限域多项式（6 个函数）

处理双变量多项式 $𝔽_p[x,y]$。

### `ocas_poly_fp_create`

**签名**：
```c
OcasPolyFp *ocas_poly_fp_create(const char *input, const char *prime, int *err);
```

**功能**：在素数域 $𝔽_p$ 上创建双变量多项式。字符串中的系数自动取模 $p$。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `input` | `const char*` | 多项式字符串 |
| `prime` | `const char*` | 素数 $p$ 的十进制字符串 |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasPolyFp*` 句柄或 `NULL`。

**内存**：调用方用 `ocas_poly_fp_free()` 释放。

**示例**：
```c
OcasPolyFp *p = ocas_poly_fp_create("x^2 + 2*x + 1", "7", NULL);
// 系数在 𝔽₇ 中：x² + 2x + 1
ocas_poly_fp_free(p);
```

---

### `ocas_poly_fp_free`

**签名**：
```c
void ocas_poly_fp_free(OcasPolyFp *poly);
```

**功能**：释放有限域多项式句柄。`NULL` 安全。

---

### `ocas_poly_fp_clone`

**签名**：
```c
OcasPolyFp *ocas_poly_fp_clone(const OcasPolyFp *poly);
```

**功能**：克隆有限域多项式。返回新句柄。

---

### `ocas_poly_fp_degree`

**签名**：
```c
size_t ocas_poly_fp_degree(const OcasPolyFp *poly);
```

**功能**：返回多项式的总次数。`NULL` 返回 `0`。

---

### `ocas_poly_fp_to_string`

**签名**：
```c
char *ocas_poly_fp_to_string(const OcasPolyFp *poly, int *err);
```

**功能**：返回多项式的堆分配字符串。调用方用 `ocas_string_free()` 释放。

---

### `ocas_poly_fp_factor`

**签名**：
```c
int ocas_poly_fp_factor(const OcasPolyFp *poly, OcasPolyFactorArray *out, int *err);
```

**功能**：对有限域多项式进行因式分解。返回值和内存管理同 `ocas_poly_z_factor`，但 `out.factors[i].poly` 需转换为 `OcasPolyFp*`。

**示例**：
```c
OcasPolyFp *p = ocas_poly_fp_create("x^2 - 1", "5", NULL);
OcasPolyFactorArray arr = {0};
if (ocas_poly_fp_factor(p, &arr, NULL) == OCAS_OK) {
    for (size_t i = 0; i < arr.len; i++) {
        OcasPolyFp *f = (OcasPolyFp *)arr.factors[i].poly;
        char *s = ocas_poly_fp_to_string(f, NULL);
        printf("因子: %s (重数 %zu)\n", s, arr.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_fp_free(f);
    }
    ocas_poly_factor_array_free(&arr);
}
ocas_poly_fp_free(p);
```

---

## 5. 因子数组（1 个函数）

### `ocas_poly_factor_array_free`

**签名**：
```c
void ocas_poly_factor_array_free(OcasPolyFactorArray *arr);
```

**功能**：释放因子数组结构本身及内部的 `OcasPolyFactor` 对象，但**不释放**各因子的多项式句柄。

> **⚠️ 内存陷阱**：必须先用 `ocas_poly_z_free()` / `ocas_poly_fp_free()` 释放每个 `arr->factors[i].poly`，再调用此函数。顺序颠倒会导致 use-after-free。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `arr` | `OcasPolyFactorArray*` | 由 `ocas_poly_z_factor` 或 `ocas_poly_fp_factor` 填充的数组 |

---

## 6. 代数数域（9 个函数）

处理代数数域 $ℚ(α)$ 上的多项式。支持 Trager 因式分解算法。

### 数据格式

**极小多项式**：单变量字符串，变量必须为 `x`，如 `"x^2 - 2"` 表示 $ℚ(\sqrt{2})$。

**系数列表**：分号分隔（常数项在前），每个系数是逗号分隔的有理数列表（$α$-多项式，升序）。例如在 $ℚ(\sqrt{2})$ 上：
- `"-2;0;1"` — $x^2 - 2$（所有系数在基域中）
- `"0;0;0,1"` — $x^2 - α$（$x^2$ 系数为 $0 + 1 \cdot α$）

单个有理数可省略逗号。

### `ocas_algebraic_field_create`

**签名**：
```c
OcasAlgebraicField *ocas_algebraic_field_create(const char *min_poly, int *err);
```

**功能**：从首一极小多项式创建代数数域。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `min_poly` | `const char*` | 极小多项式字符串，如 `"x^2 - 2"` |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasAlgebraicField*` 句柄或 `NULL`。

**内存**：调用方用 `ocas_algebraic_field_free()` 释放。

**示例**：
```c
OcasAlgebraicField *q_sqrt2 = ocas_algebraic_field_create("x^2 - 2", NULL);
printf("扩张次数: %zu\n", ocas_algebraic_field_degree(q_sqrt2));
// 输出：扩张次数: 2
ocas_algebraic_field_free(q_sqrt2);
```

---

### `ocas_algebraic_field_free`

**签名**：
```c
void ocas_algebraic_field_free(OcasAlgebraicField *field);
```

**功能**：释放代数数域句柄。`NULL` 安全。

---

### `ocas_algebraic_field_degree`

**签名**：
```c
size_t ocas_algebraic_field_degree(const OcasAlgebraicField *field);
```

**功能**：返回扩张次数 $\deg(m) = [ℚ(α):ℚ]$。`NULL` 返回 `0`。

---

### `ocas_algebraic_poly_create`

**签名**：
```c
OcasAlgebraicPoly *ocas_algebraic_poly_create(const OcasAlgebraicField *field,
                                              const char *coeffs, int *err);
```

**功能**：在代数数域上创建多项式。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `field` | `const OcasAlgebraicField*` | 非空数域句柄 |
| `coeffs` | `const char*` | 系数列表字符串（见上文格式） |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasAlgebraicPoly*` 句柄或 `NULL`。

**示例**：
```c
OcasAlgebraicField *fld = ocas_algebraic_field_create("x^2 - 2", NULL);
// 创建多项式 x² - α（即 x² - √2）
OcasAlgebraicPoly *p = ocas_algebraic_poly_create(fld, "0;0;0,1", NULL);
printf("次数: %zu\n", ocas_algebraic_poly_degree(p));
ocas_algebraic_poly_free(p);
ocas_algebraic_field_free(fld);
```

---

### `ocas_algebraic_poly_free`

**签名**：
```c
void ocas_algebraic_poly_free(OcasAlgebraicPoly *poly);
```

**功能**：释放代数域多项式句柄。`NULL` 安全。

---

### `ocas_algebraic_poly_degree`

**签名**：
```c
size_t ocas_algebraic_poly_degree(const OcasAlgebraicPoly *poly);
```

**功能**：返回多项式次数。零多项式返回 `0`。

---

### `ocas_algebraic_poly_to_string`

**签名**：
```c
char *ocas_algebraic_poly_to_string(const OcasAlgebraicPoly *poly, int *err);
```

**功能**：返回多项式的堆分配字符串。格式为 `[c0] + [c1]*x + [c2]*x^2 + ...`，其中 `[ci]` 是逗号分隔的 $α$-多项式有理数。调用方用 `ocas_string_free()` 释放。

---

### `ocas_algebraic_poly_factor`

**签名**：
```c
int ocas_algebraic_poly_factor(const OcasAlgebraicPoly *poly,
                               OcasAlgebraicFactorArray *out, int *err);
```

**功能**：使用 Trager 算法在代数数域上因式分解多项式。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `poly` | `const OcasAlgebraicPoly*` | 非空多项式句柄 |
| `out` | `OcasAlgebraicFactorArray*` | 输出数组（成功时填充） |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OCAS_OK`（0）成功。

**内存**：`out->factors[i].poly`（`void*`）需转换为 `OcasAlgebraicPoly*` 并用 `ocas_algebraic_poly_free()` 释放。数组用 `ocas_algebraic_factor_array_free()` 释放。

**示例**：
```c
OcasAlgebraicField *fld = ocas_algebraic_field_create("x^2 - 2", NULL);
OcasAlgebraicPoly *p = ocas_algebraic_poly_create(fld, "-2;0;1", NULL);
OcasAlgebraicFactorArray arr = {0};
if (ocas_algebraic_poly_factor(p, &arr, NULL) == OCAS_OK) {
    for (size_t i = 0; i < arr.len; i++) {
        OcasAlgebraicPoly *f = (OcasAlgebraicPoly *)arr.factors[i].poly;
        char *s = ocas_algebraic_poly_to_string(f, NULL);
        printf("因子: %s (重数 %zu)\n", s, arr.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_algebraic_poly_free(f);
    }
    ocas_algebraic_factor_array_free(&arr);
}
ocas_algebraic_poly_free(p);
ocas_algebraic_field_free(fld);
```

---

### `ocas_algebraic_factor_array_free`

**签名**：
```c
void ocas_algebraic_factor_array_free(OcasAlgebraicFactorArray *arr);
```

**功能**：释放因子数组的存储。**不释放**各因子的多项式句柄——必须先单独释放。

---

## 7. 数值积分（6 个函数）

Vegas 自适应蒙特卡洛积分器。

### 积分子函数签名

```c
typedef double (*ocas_integrand_t)(double x, void *user_data);
```

`user_data` 从调用方原样传递到积分子。

> **⚠️ 注意**：当前 C API 的积分子只接收采样点的**第一个坐标**（一维积分友好接口）。即使以 `n_dims > 1` 创建积分器，积分子的 `x` 参数也仅代表第一个维度；完整的多维 `const double*`/`size_t` 签名留待后续版本。

### `ocas_vegas_create`

**签名**：
```c
OcasVegas *ocas_vegas_create(size_t n_dims, const OcasVegasOptions *opts, int *err);
```

**功能**：创建 `n_dims` 维 Vegas 积分器。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `n_dims` | `size_t` | 积分维度 |
| `opts` | `const OcasVegasOptions*` | 可为 `NULL`（使用默认值） |
| `err` | `int*` | 可为 `NULL` |

**`OcasVegasOptions` 字段**：
| 字段 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `n_bins` | `size_t` | 64 | 每维度 bin 数 |
| `n_samples` | `size_t` | 10000 | 每次迭代采样数 |
| `iterations` | `size_t` | 10 | 自适应迭代次数 |
| `learning_rate` | `double` | 1.5 | 网格平滑/学习率 |
| `seed` | `uint64_t` | 0x0C45 | 随机数种子 |

**返回值**：`OcasVegas*` 句柄或 `NULL`。

**内存**：调用方用 `ocas_vegas_free()` 释放。

**示例**：
```c
OcasVegasOptions opts = {
    .n_bins = 50, .n_samples = 20000, .iterations = 8,
    .learning_rate = 1.5, .seed = 42
};
int err = 0;
OcasVegas *v = ocas_vegas_create(1, &opts, &err);
```

---

### `ocas_vegas_free`

**签名**：
```c
void ocas_vegas_free(OcasVegas *v);
```

**功能**：释放 Vegas 积分器句柄。`NULL` 安全。

---

### `ocas_vegas_integrate`

**签名**：
```c
OcasIntegrateResult ocas_vegas_integrate(OcasVegas *v, ocas_integrand_t f,
                                         void *user_data, int *err);
```

**功能**：在单位超立方体上对 `f` 积分（积分子仅接收第一个坐标，见上文注意）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `v` | `OcasVegas*` | 非空积分器句柄 |
| `f` | `ocas_integrand_t` | 积分子函数 |
| `user_data` | `void*` | 传递给积分子的用户数据 |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasIntegrateResult { double integral; double error; }`。

**示例**：
```c
static double f(double x, void *ud) { return x * x; }

OcasIntegrateResult r = ocas_vegas_integrate(v, f, NULL, &err);
printf("∫₀¹ x² dx ≈ %g ± %g\n", r.integral, r.error);
// 输出（示例值，随种子与采样数而定）：∫₀¹ x² dx ≈ 0.333 ± 0.001
```

---

### `ocas_vegas_result`

**签名**：
```c
OcasIntegrateResult ocas_vegas_result(const OcasVegas *v);
```

**功能**：返回最近一次 `ocas_vegas_integrate` 的累计估计值和误差。

---

### `ocas_vegas_iterations`

**签名**：
```c
size_t ocas_vegas_iterations(const OcasVegas *v);
```

**功能**：返回已完成的迭代次数。`NULL` 返回 `0`。

---

### `ocas_integrate_1d`

**签名**：
```c
OcasIntegrateResult ocas_integrate_1d(ocas_integrand_t f, void *user_data,
                                      double a, double b,
                                      const OcasVegasOptions *opts, int *err);
```

**功能**：一维数值积分的便捷函数——在 $[a, b]$ 上使用 Vegas 积分。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `f` | `ocas_integrand_t` | 积分子 |
| `user_data` | `void*` | 用户数据 |
| `a` | `double` | 积分下限 |
| `b` | `double` | 积分上限 |
| `opts` | `const OcasVegasOptions*` | 可为 `NULL` |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasIntegrateResult`。

**示例**：
```c
OcasIntegrateResult r = ocas_integrate_1d(f, NULL, 0.0, 1.0, NULL, &err);
printf("结果: %g ± %g\n", r.integral, r.error);
```

---

## 8. ODE（3 个函数）

符号常微分方程求解。所有输入输出均为字符串。

### `ocas_ode_classify`

**签名**：
```c
char *ocas_ode_classify(const char *equation, const char *func,
                        const char *var, int *err);
```

**功能**：分类 ODE 并返回适用的求解方法名列表（逗号分隔）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `equation` | `const char*` | ODE 等于零的表达式，如 `"Derivative(y(x), x) - y(x)"` |
| `func` | `const char*` | 未知函数名，如 `"y"` |
| `var` | `const char*` | 自变量，如 `"x"` |
| `err` | `int*` | 可为 `NULL` |

**返回值**：堆分配字符串如 `"LinearFirst,PowerSeries"`，调用方用 `ocas_string_free()` 释放。失败返回 `NULL`。

**示例**：
```c
char *methods = ocas_ode_classify("Derivative(y(x),x) - y(x)", "y", "x", NULL);
printf("可用方法: %s\n", methods);
// 输出：可用方法: LinearFirst,PowerSeries
ocas_string_free(methods);
```

---

### `ocas_ode_dsolve`

**签名**：
```c
char *ocas_ode_dsolve(const char *equation, const char *func,
                      const char *var, const char *hint, int *err);
```

**功能**：符号求解 ODE。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `equation` | `const char*` | ODE 表达式 |
| `func` | `const char*` | 未知函数名 |
| `var` | `const char*` | 自变量 |
| `hint` | `const char*` | 可为 `NULL`（自动分类）；或 `ocas_ode_classify` 返回的方法名之一 |
| `err` | `int*` | 可为 `NULL` |

**返回值**：解字符串如 `"y = C1*exp(x)"` 或 `"unsolved"`。调用方用 `ocas_string_free()` 释放。

**示例**：
```c
char *sol = ocas_ode_dsolve("Derivative(y(x),x) - y(x)", "y", "x", NULL, NULL);
printf("解: %s\n", sol);  // 输出：解: y = C1*exp(x)
ocas_string_free(sol);
```

---

### `ocas_ode_dsolve_ivp`

**签名**：
```c
char *ocas_ode_dsolve_ivp(const char *equation, const char *func,
                          const char *var, const char *y0,
                          const char *y1, int *err);
```

**功能**：通过 Laplace 变换求解一阶或二阶线性常系数 ODE 的初值问题。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `equation` | `const char*` | ODE 表达式 |
| `func` | `const char*` | 未知函数名 |
| `var` | `const char*` | 自变量 |
| `y0` | `const char*` | $y(0)$ 的表达式字符串，如 `"1"` |
| `y1` | `const char*` | $y'(0)$（可为 `NULL`，仅二阶需要） |
| `err` | `int*` | 可为 `NULL` |

**返回值**：无自由常数的显式解字符串。调用方用 `ocas_string_free()` 释放。

**示例**：
```c
// y'' + y = 0, y(0) = 1, y'(0) = 0 → y = cos(x)
char *sol = ocas_ode_dsolve_ivp(
    "Derivative(y(x),x,2) + y(x)", "y", "x", "1", "0", NULL);
printf("IVP 解: %s\n", sol);
ocas_string_free(sol);
```

---

## 9. 数论（11 个函数）

任意精度整数通过十进制字符串传递。字符串结果由调用方用 `ocas_string_free()` 释放。

### `ocas_ntheory_factorint`

**签名**：
```c
char *ocas_ntheory_factorint(const char *n, int *err_out);
```

**功能**：对 $|n|$ 进行素因数分解。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `const char*` | 十进制整数字符串 |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：`"p1:e1,p2:e2,..."`（升序），负数开头含 `"-1:1"`。失败返回 `NULL`。

**示例**：
```c
char *f = ocas_ntheory_factorint("360", NULL);
printf("360 = %s\n", f);  // 输出：360 = 2:3,3:2,5:1
ocas_string_free(f);

f = ocas_ntheory_factorint("-12", NULL);
printf("-12 = %s\n", f);  // 输出：-12 = -1:1,2:2,3:1
ocas_string_free(f);
```

---

### `ocas_ntheory_isprime`

**签名**：
```c
int ocas_ntheory_isprime(const char *n, int *err_out);
```

**功能**：BPSW 概率素性测试。

**返回值**：`1` = (很可能) 素数，`0` = 合数，`-1` = 错误。

**示例**：
```c
int r = ocas_ntheory_isprime("997", NULL);
printf("997 是素数？%s\n", r == 1 ? "是" : "否");
// 输出：997 是素数？是
```

---

### `ocas_ntheory_nextprime`

**签名**：
```c
char *ocas_ntheory_nextprime(const char *n, int *err_out);
```

**功能**：返回严格大于 `n` 的最小素数。

**返回值**：十进制字符串或 `NULL`。

**示例**：
```c
char *p = ocas_ntheory_nextprime("100", NULL);
printf("100 之后的素数: %s\n", p);  // 输出：101
ocas_string_free(p);
```

---

### `ocas_ntheory_discrete_log`

**签名**：
```c
char *ocas_ntheory_discrete_log(const char *p, const char *base,
                                const char *target, int *err_out);
```

**功能**：求解 $\text{base}^x \equiv \text{target} \pmod{p}$。素数 $p$ 使用 Pohlig–Hellman，否则使用 BSGS。

**返回值**：对数的十进制字符串，或 `NULL`（无解时设 `OCAS_ERROR_RUNTIME`）。

**示例**：
```c
char *x = ocas_ntheory_discrete_log("7", "3", "5", NULL);
// 3^x ≡ 5 (mod 7) → x = 5
printf("离散对数: %s\n", x);
ocas_string_free(x);
```

---

### `ocas_ntheory_crt`

**签名**：
```c
char *ocas_ntheory_crt(const char *moduli, const char *residues, int *err_out);
```

**功能**：中国剩余定理。模数无需互素；不一致返回 `NULL`。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `moduli` | `const char*` | 逗号分隔模数，如 `"3,5,7"` |
| `residues` | `const char*` | 逗号分隔余数，如 `"2,3,2"` |

**返回值**：`"r,m"`（$r \equiv \text{residues}[i] \pmod{\text{moduli}[i]}$），或 `NULL`。

**示例**：
```c
char *r = ocas_ntheory_crt("3,5,7", "2,3,2", NULL);
printf("CRT: %s\n", r);  // 输出：CRT: 23,105
// 23 ≡ 2 (mod 3), 23 ≡ 3 (mod 5), 23 ≡ 2 (mod 7)
ocas_string_free(r);
```

---

### `ocas_ntheory_jacobi`

**签名**：
```c
int ocas_ntheory_jacobi(const char *a, const char *n, int *err_out);
```

**功能**：计算 Jacobi 符号 $(a/n)$，要求 $n$ 为正奇数。

**返回值**：$-1$、$0$ 或 $1$；$-2$ 表示输入无效。

---

### `ocas_ntheory_totient`

**签名**：
```c
char *ocas_ntheory_totient(const char *n, int *err_out);
```

**功能**：Euler totient 函数 $\varphi(n)$。

**返回值**：十进制字符串或 `NULL`。

**示例**：
```c
char *t = ocas_ntheory_totient("12", NULL);
printf("φ(12) = %s\n", t);  // 输出：φ(12) = 4
ocas_string_free(t);
```

---

### `ocas_ntheory_mobius`

**签名**：
```c
int ocas_ntheory_mobius(const char *n, int *err_out);
```

**功能**：Möbius 函数 $\mu(n)$。

**返回值**：$-1$、$0$ 或 $1$；$-2$ 表示错误。

---

### `ocas_ntheory_divisor_count`

**签名**：
```c
char *ocas_ntheory_divisor_count(const char *n, int *err_out);
```

**功能**：正因数个数 $\tau(n)$。

**返回值**：十进制字符串或 `NULL`。

**示例**：
```c
char *d = ocas_ntheory_divisor_count("12", NULL);
printf("τ(12) = %s\n", d);  // 输出：τ(12) = 6
// 12 的因数：1,2,3,4,6,12
ocas_string_free(d);
```

---

### `ocas_ntheory_divisor_sigma`

**签名**：
```c
char *ocas_ntheory_divisor_sigma(const char *n, uint32_t k, int *err_out);
```

**功能**：正因数的 $k$ 次幂和 $\sigma_k(n)$。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `n` | `const char*` | 正整数 |
| `k` | `uint32_t` | 幂次 |
| `err_out` | `int*` | 可为 `NULL` |

**返回值**：十进制字符串或 `NULL`。

**示例**：
```c
char *s = ocas_ntheory_divisor_sigma("6", 1, NULL);
printf("σ₁(6) = %s\n", s);  // 输出：σ₁(6) = 12
// 1+2+3+6 = 12
ocas_string_free(s);
```

---

### `ocas_ntheory_liouville`

**签名**：
```c
int ocas_ntheory_liouville(const char *n, int *err_out);
```

**功能**：Liouville 函数 $\lambda(n) = (-1)^{\Omega(n)}$。

**返回值**：$-1$、$0$ 或 $1$；$-2$ 表示错误。

---

### C++ 包装：`ocas::ntheory::*`

所有数论函数均有 C++ 内联包装，自动管理字符串内存：

```cpp
#include <ocas.h>
#include <ocas.hpp>

std::string f = ocas::ntheory::factorint("360");
// f == "2:3,3:2,5:1"

bool prime = ocas::ntheory::isprime("997");
// prime == true

std::string p = ocas::ntheory::nextprime("100");
// p == "101"

std::string tot = ocas::ntheory::totient("12");
// tot == "4"
```

| C++ 函数 | C 函数 | 说明 |
|---|---|---|
| `factorint(n)` | `ocas_ntheory_factorint` | 素因数分解 |
| `isprime(n)` | `ocas_ntheory_isprime` | 返回 `bool`，错误抛异常 |
| `nextprime(n)` | `ocas_ntheory_nextprime` | 下一个素数 |
| `discrete_log(p, base, target)` | `ocas_ntheory_discrete_log` | 离散对数，无解抛异常 |
| `crt(moduli, residues)` | `ocas_ntheory_crt` | 中国剩余定理 |
| `jacobi(a, n)` | `ocas_ntheory_jacobi` | Jacobi 符号，无效输入抛异常 |
| `totient(n)` | `ocas_ntheory_totient` | Euler totient |
| `mobius(n)` | `ocas_ntheory_mobius` | Möbius 函数，无效输入抛异常 |
| `divisor_count(n)` | `ocas_ntheory_divisor_count` | 因数个数 |
| `divisor_sigma(n, k=1)` | `ocas_ntheory_divisor_sigma` | 因数幂和 |
| `liouville(n)` | `ocas_ntheory_liouville` | Liouville 函数，无效输入抛异常 |

所有包装在 C API 返回错误时抛出 `ocas::Error` 异常。

---

## 10. 张量（12 个函数）

命名指标张量，支持指标缩并和对称性。

### 槽位字符串格式

槽位以分号分隔，每个条目为 `label,position`，其中 `position` 为 `"upper"` 或 `"lower"`（别名：`up`/`down`/`contravariant`/`covariant`）。

例如：`"i,upper;j,lower"` 表示两个槽位 $i^j{}_{}$。

### `ocas_tensor_create`

**签名**：
```c
OcasTensor *ocas_tensor_create(const char *name, const char *slots,
                               const char *symmetry, int *err);
```

**功能**：创建张量。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `name` | `const char*` | 张量名，如 `"T"` |
| `slots` | `const char*` | 槽位字符串 |
| `symmetry` | `const char*` | 可为 `NULL`（= `"none"`）；或 `"none"`/`"symmetric"`/`"antisymmetric"` |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasTensor*` 句柄或 `NULL`。

**示例**：
```c
OcasTensor *T = ocas_tensor_create("T", "i,upper;j,lower", "symmetric", NULL);
printf("秩: %zu\n", ocas_tensor_rank(T));  // 输出：秩: 2
ocas_tensor_free(T);
```

---

### `ocas_tensor_free`

**签名**：
```c
void ocas_tensor_free(OcasTensor *t);
```

**功能**：释放张量句柄。`NULL` 安全。

---

### `ocas_tensor_name`

**签名**：
```c
char *ocas_tensor_name(const OcasTensor *t, int *err);
```

**功能**：返回张量名称的堆分配字符串。调用方用 `ocas_string_free()` 释放。

---

### `ocas_tensor_rank`

**签名**：
```c
size_t ocas_tensor_rank(const OcasTensor *t);
```

**功能**：返回张量的秩（槽位数）。`NULL` 返回 `0`。

---

### `ocas_tensor_symmetry`

**签名**：
```c
int ocas_tensor_symmetry(const OcasTensor *t);
```

**功能**：返回对称性代码。

**返回值**：`0` = none，`1` = symmetric，`2` = antisymmetric，`-1` = null 句柄。

---

### `ocas_tensor_to_string`

**签名**：
```c
char *ocas_tensor_to_string(const OcasTensor *t, int *err);
```

**功能**：返回张量的字符串表示 `name(slot, slot, ...)`。调用方用 `ocas_string_free()` 释放。

---

### `ocas_tensor_symmetrise_sign`

**签名**：
```c
int64_t ocas_tensor_symmetrise_sign(const OcasTensor *t);
```

**功能**：返回反对称化符号（$+1$ 或 $-1$）。`NULL` 返回 `0`。

---

### `ocas_tensor_contract`

**签名**：
```c
int ocas_tensor_contract(const OcasTensor *a, const OcasTensor *b,
                         OcasTensorContraction *out, int *err);
```

**功能**：对两个张量的共享哑指标（相同标签、相反变异性）求和缩并。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `a` | `const OcasTensor*` | 非空张量 |
| `b` | `const OcasTensor*` | 非空张量 |
| `out` | `OcasTensorContraction*` | 输出结构体 |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OCAS_OK`（0）成功。

**`OcasTensorContraction` 字段**：
| 字段 | 类型 | 说明 |
|---|---|---|
| `kind` | `int` | `0` = product（有自由指标），`1` = scalar（完全缩并） |
| `tensors` | `OcasTensor**` | `kind == 0` 时有效，张量句柄数组 |
| `n_tensors` | `size_t` | `tensors` 数组长度 |
| `scalar_str` | `char*` | `kind == 1` 时有效，标量结果字符串 |

**内存**：
- `kind == 0`：每个 `tensors[i]` 用 `ocas_tensor_free()` 释放，数组用 `ocas_tensor_contraction_free()` 释放。
- `kind == 1`：`scalar_str` 用 `ocas_string_free()` 释放，结构体用 `ocas_tensor_contraction_free()` 释放。

**示例**：
```c
OcasTensor *A = ocas_tensor_create("A", "i,upper;j,lower", NULL, NULL);
OcasTensor *B = ocas_tensor_create("B", "j,upper;k,lower", NULL, NULL);
OcasTensorContraction result = {0};
if (ocas_tensor_contract(A, B, &result, NULL) == OCAS_OK) {
    if (result.kind == 0) {
        // 自由指标剩余
        for (size_t i = 0; i < result.n_tensors; i++) {
            char *s = ocas_tensor_to_string(result.tensors[i], NULL);
            printf("结果张量: %s\n", s);
            ocas_string_free(s);
            ocas_tensor_free(result.tensors[i]);
        }
    }
    ocas_tensor_contraction_free(&result);
}
ocas_tensor_free(B);
ocas_tensor_free(A);
```

---

### `ocas_tensor_contraction_free`

**签名**：
```c
void ocas_tensor_contraction_free(OcasTensorContraction *c);
```

**功能**：释放缩并结果的 `tensors` 数组和 `scalar_str`。**不释放**各张量句柄。`NULL` 安全。

---

### `ocas_tensor_canonicalize`

**签名**：
```c
char *ocas_tensor_canonicalize(const char *expr_str, const char *specs_str,
                                const char *groups_str, int *err);
```

**功能**：通过图同构规范化张量表达式（0.22.0）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `expr_str` | `const char*` | 张量积字符串，如 `"T(i,j)*U(j,k)"` |
| `specs_str` | `const char*` | 逗号分隔的 `name:sym` 对，如 `"T:none,U:none"` |
| `groups_str` | `const char*` | 可为 `NULL`；或逗号分隔的 `label:group` 对 |
| `err` | `int*` | 可为 `NULL` |

**返回值**：堆分配的规范形式字符串。调用方用 `ocas_string_free()` 释放。

**示例**：
```c
char *canon = ocas_tensor_canonicalize(
    "T(i,j)*T(j,i)", "T:symmetric", NULL, NULL);
printf("规范形式: %s\n", canon);
ocas_string_free(canon);
```

---

### `ocas_young_project`

**签名**：
```c
char *ocas_young_project(const char *expr_str, const char *tableau_str, int *err);
```

**功能**：对张量表达式应用 Young 投影子（0.22.0）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `expr_str` | `const char*` | 张量表达式，如 `"f(a,b)"` |
| `tableau_str` | `const char*` | 逗号分隔的行长度，如 `"2,1"` 表示 □□/□ |
| `err` | `int*` | 可为 `NULL` |

**返回值**：堆分配的展开表达式字符串。调用方用 `ocas_string_free()` 释放。

**示例**：
```c
// 全反对称投影（(1,1,1) Young 表，三指标）
char *proj = ocas_young_project("f(a,b,c)", "1,1,1", NULL);
printf("投影结果: %s\n", proj);
ocas_string_free(proj);
```

---

### `ocas_tensor_refresh_dummies`

**签名**：
```c
char *ocas_tensor_refresh_dummies(const char *expr_str, const char *specs_str, int *err);
```

**功能**：重命名张量表达式中的哑指标（恰好出现两次的标签）为 `d0`、`d1`、… 以避免冲突（0.22.0）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `expr_str` | `const char*` | 张量表达式 |
| `specs_str` | `const char*` | `name:sym` 对 |
| `err` | `int*` | 可为 `NULL` |

**返回值**：堆分配字符串。调用方用 `ocas_string_free()` 释放。

---

### C++ 包装：`ocas::tensor::*`

```cpp
std::string canon = ocas::tensor::canonicalize("T(i,j)*T(j,i)", "T:symmetric");
std::string proj = ocas::tensor::young_project("f(a,b,c)", "2,1");
std::string fresh = ocas::tensor::refresh_dummies("T(i,j)*U(j,i)", "T:none,U:none");
```

所有包装在错误时抛出 `ocas::Error`。

---

## 11. 双数（14 个函数）

超对偶数（前向自动微分），仅支持 `Rational` 系数。多项式/有理算术。

### 系数字符串格式

有理数以 `"num"`（分母为 1）或 `"num/den"` 传递。返回字符串使用相同格式。

### `ocas_dual_shape_new`

**签名**：
```c
OcasDualShape *ocas_dual_shape_new(size_t n_vars, int *err);
```

**功能**：创建一阶形状，跟踪 `n_vars` 个变量的导数。

**返回值**：`OcasDualShape*` 句柄或 `NULL`。

**示例**：
```c
OcasDualShape *shape = ocas_dual_shape_new(2, NULL);  // 2 个变量
printf("变量数: %zu, 分量数: %zu\n",
       ocas_dual_shape_n_vars(shape),
       ocas_dual_shape_n_components(shape));
// 输出：变量数: 2, 分量数: 3
ocas_dual_shape_free(shape);
```

---

### `ocas_dual_shape_free`

**签名**：
```c
void ocas_dual_shape_free(OcasDualShape *s);
```

**功能**：释放形状句柄。`NULL` 安全。

---

### `ocas_dual_shape_n_vars`

**签名**：
```c
size_t ocas_dual_shape_n_vars(const OcasDualShape *s);
```

**功能**：返回微分变量数。`NULL` 返回 `0`。

---

### `ocas_dual_shape_n_components`

**签名**：
```c
size_t ocas_dual_shape_n_components(const OcasDualShape *s);
```

**功能**：返回总分量数（值 + 各偏导）。`NULL` 返回 `0`。

---

### `ocas_dual_variable`

**签名**：
```c
OcasHyperDual *ocas_dual_variable(const OcasDualShape *shape, size_t i,
                                  const char *coeff, int *err);
```

**功能**：创建独立变量 $x_i = \text{coeff}$（对变量 $i$ 的导数为 1）。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `shape` | `const OcasDualShape*` | 非空形状句柄 |
| `i` | `size_t` | 变量索引（$0 \le i < \text{n\_vars}$） |
| `coeff` | `const char*` | 系数字符串 `"num"` 或 `"num/den"` |
| `err` | `int*` | 可为 `NULL` |

**返回值**：`OcasHyperDual*` 句柄或 `NULL`。

---

### `ocas_dual_constant`

**签名**：
```c
OcasHyperDual *ocas_dual_constant(const OcasDualShape *shape,
                                  const char *coeff, int *err);
```

**功能**：创建常量双数（所有导数为零）。

---

### `ocas_hyperdual_free`

**签名**：
```c
void ocas_hyperdual_free(OcasHyperDual *d);
```

**功能**：释放超对偶数句柄。`NULL` 安全。

---

### `ocas_dual_value`

**签名**：
```c
char *ocas_dual_value(const OcasHyperDual *d, int *err);
```

**功能**：返回标量值分量的堆分配字符串。调用方用 `ocas_string_free()` 释放。

---

### `ocas_dual_deriv`

**签名**：
```c
char *ocas_dual_deriv(const OcasHyperDual *d, size_t i, int *err);
```

**功能**：返回关于变量 $i$ 的偏导数的堆分配字符串。形状无 $i$ 的一阶分量时返回 `NULL`。

**示例**：
```c
OcasDualShape *shape = ocas_dual_shape_new(2, NULL);
OcasHyperDual *x = ocas_dual_variable(shape, 0, "3", NULL);   // x₀ = 3
OcasHyperDual *y = ocas_dual_variable(shape, 1, "5", NULL);   // x₁ = 5
OcasHyperDual *prod = ocas_dual_mul(x, y, NULL);               // f = x₀·x₁

char *val = ocas_dual_value(prod, NULL);
char *df_dx = ocas_dual_deriv(prod, 0, NULL);  // ∂f/∂x₀ = x₁ = 5
char *df_dy = ocas_dual_deriv(prod, 1, NULL);  // ∂f/∂x₁ = x₀ = 3

printf("f = %s, ∂f/∂x₀ = %s, ∂f/∂x₁ = %s\n", val, df_dx, df_dy);
// 输出：f = 15, ∂f/∂x₀ = 5, ∂f/∂x₁ = 3

ocas_string_free(val);
ocas_string_free(df_dx);
ocas_string_free(df_dy);
ocas_hyperdual_free(prod);
ocas_hyperdual_free(y);
ocas_hyperdual_free(x);
ocas_dual_shape_free(shape);
```

---

### `ocas_dual_add`

**签名**：
```c
OcasHyperDual *ocas_dual_add(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**功能**：计算 $a + b$。两个操作数必须共享同一形状。

**返回值**：新句柄或 `NULL`。

---

### `ocas_dual_sub`

**签名**：
```c
OcasHyperDual *ocas_dual_sub(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**功能**：计算 $a - b$。

---

### `ocas_dual_mul`

**签名**：
```c
OcasHyperDual *ocas_dual_mul(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**功能**：计算 $a \times b$。通过乘积法则自动传播导数。

---

### `ocas_dual_div`

**签名**：
```c
OcasHyperDual *ocas_dual_div(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**功能**：计算 $a / b$。若 $b$ 的值分量为零，设置 `OCAS_ERROR_DIVISION_BY_ZERO` 并返回 `NULL`。

---

### `ocas_dual_neg`

**签名**：
```c
OcasHyperDual *ocas_dual_neg(const OcasHyperDual *a, int *err);
```

**功能**：计算 $-a$。

---

## 12. Gröbner（8 个函数）

Gröbner 基计算和多项式系统求解。

> **⚠️ 非标准约定**
>
> - 这些函数**未包含在 `ocas.h` 头文件中**，需自行声明原型。
> - 错误码使用 `-1` 而非 `OCAS_*` 常量。
> - `OcasGroebnerBasis` 和 `OcasSystemSolution` 句柄**未标记** `#[repr(C)]`。

### 多项式数据格式

多项式通过系数数组传递：每个多项式由以下部分指定：
- `n_vars` — 变量数
- `n_terms` — 项数
- `exponents` — 展平的指数矩阵（`n_terms × n_vars`，行优先）
- `coeff_nums` / `coeff_dens` — 有理系数的分子/分母数组

### `ocas_groebner_basis`

**签名**：
```c
OcasGroebnerBasis *ocas_groebner_basis(
    size_t n_polys,
    const size_t *n_vars_array,
    const size_t *n_terms_array,
    const size_t *exponents,
    const int64_t *coeff_nums,
    const int64_t *coeff_dens,
    int32_t algorithm,
    int32_t *err
);
```

**功能**：从多项式数据数组计算 Gröbner 基。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `n_polys` | `size_t` | 多项式个数 |
| `n_vars_array` | `const size_t*` | 每个多项式的变量数 |
| `n_terms_array` | `const size_t*` | 每个多项式的项数 |
| `exponents` | `const size_t*` | 展平指数矩阵 |
| `coeff_nums` | `const int64_t*` | 有理系数分子 |
| `coeff_dens` | `const int64_t*` | 有理系数分母 |
| `algorithm` | `int32_t` | 算法选择：`0` = Auto, `1` = F4, `2` = F5, `3` = Buchberger |
| `err` | `int32_t*` | 错误码（`-1` 表示失败）。**注意**：与其余 C API 不同，此指针不可为 `NULL`——为 `NULL` 时函数直接返回 `NULL` 且不计算 |

**返回值**：`OcasGroebnerBasis*` 句柄或 `NULL`。

**内存**：调用方用 `ocas_groebner_basis_free()` 释放。

**示例**：
```c
// 计算 {x² + y - 1, x + y² - 1} 的 Gröbner 基
size_t n_vars[] = {2, 2};
size_t n_terms[] = {2, 2};
// x² + y - 1: exponents = [[2,0],[0,1],[0,0]]
// x + y² - 1: exponents = [[1,0],[0,2],[0,0]]
size_t exps[] = {2,0, 0,1, 0,0,
                 1,0, 0,2, 0,0};
int64_t nums[] = {1, 1, -1,  1, 1, -1};
int64_t dens[] = {1, 1,  1,  1, 1,  1};

int32_t err = 0;
OcasGroebnerBasis *gb = ocas_groebner_basis(
    2, n_vars, n_terms, exps, nums, dens, 0, &err);
if (gb) {
    printf("基有 %zu 个元素\n", ocas_groebner_basis_len(gb));
    printf("零维？%s\n", ocas_is_zero_dimensional(gb) ? "是" : "否");
    ocas_groebner_basis_free(gb);
}
```

---

### `ocas_groebner_basis_free`

**签名**：
```c
void ocas_groebner_basis_free(OcasGroebnerBasis *gb);
```

**功能**：释放 Gröbner 基句柄。`NULL` 安全。

---

### `ocas_groebner_basis_len`

**签名**：
```c
size_t ocas_groebner_basis_len(const OcasGroebnerBasis *gb);
```

**功能**：返回 Gröbner 基中的元素数。`NULL` 返回 `0`。

---

### `ocas_is_zero_dimensional`

**签名**：
```c
bool ocas_is_zero_dimensional(const OcasGroebnerBasis *gb);
```

**功能**：检查理想是否为零维。`NULL` 返回 `false`。

---

### `ocas_solve_polynomial_system`

**签名**：
```c
OcasSystemSolution *ocas_solve_polynomial_system(
    size_t n_polys,
    const size_t *n_vars_array,
    const size_t *n_terms_array,
    const size_t *exponents,
    const int64_t *coeff_nums,
    const int64_t *coeff_dens,
    int32_t algorithm,
    int32_t *err
);
```

**功能**：求解多项式方程组。返回句柄或 `NULL`。

**参数**：同 `ocas_groebner_basis`。

**内存**：调用方用 `ocas_system_solution_free()` 释放。

---

### `ocas_system_solution_count`

**签名**：
```c
size_t ocas_system_solution_count(const OcasSystemSolution *sol);
```

**功能**：返回解的个数。非零维或空解返回 `0`。`NULL` 返回 `0`。

---

### `ocas_system_solution_value`

**签名**：
```c
double ocas_system_solution_value(const OcasSystemSolution *sol,
                                  size_t sol_idx, size_t var_idx);
```

**功能**：获取特定解的特定变量值。

**参数**：
| 参数 | 类型 | 说明 |
|---|---|---|
| `sol` | `const OcasSystemSolution*` | 非空解句柄 |
| `sol_idx` | `size_t` | 解的索引 |
| `var_idx` | `size_t` | 变量的索引 |

**返回值**：`f64` 值。越界或错误返回 `0.0`。

**示例**：
```c
OcasSystemSolution *sol = ocas_solve_polynomial_system(
    2, n_vars, n_terms, exps, nums, dens, 0, &err);
if (sol) {
    size_t count = ocas_system_solution_count(sol);
    printf("找到 %zu 个解\n", count);
    for (size_t i = 0; i < count; i++) {
        printf("解 %zu: x = %g, y = %g\n", i,
               ocas_system_solution_value(sol, i, 0),
               ocas_system_solution_value(sol, i, 1));
    }
    ocas_system_solution_free(sol);
}
```

---

### `ocas_system_solution_free`

**签名**：
```c
void ocas_system_solution_free(OcasSystemSolution *sol);
```

**功能**：释放解句柄。`NULL` 安全。

---

## 完整示例

### C：符号微分与求值

```c
#include <ocas.h>
#include <stdio.h>
#include <assert.h>

int main(void) {
    int err = 0;

    // 解析表达式
    OcasExpr *f = ocas_expr_parse("x^3 - 6*x^2 + 11*x - 6", &err);
    assert(f != NULL);

    // 求导
    OcasExpr *df = ocas_expr_diff(f, "x", &err);
    assert(df != NULL);

    // 化简
    OcasExpr *simplified = ocas_expr_simplify(df, &err);
    assert(simplified != NULL);

    char *s = ocas_expr_to_string(simplified, &err);
    printf("f'(x) = %s\n", s);  // 输出：f'(x) = 11 + (-12*x) + (3*(x^2))

    // 清理
    ocas_string_free(s);
    ocas_expr_free(simplified);
    ocas_expr_free(df);
    ocas_expr_free(f);

    return 0;
}
```

编译：
```bash
gcc example.c -locas -o example
```

### C：自动微分（双数）

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // f(x,y) = x² + 2xy + y²，在 (3, 5) 处求值与偏导
    OcasDualShape *shape = ocas_dual_shape_new(2, &err);
    OcasHyperDual *x = ocas_dual_variable(shape, 0, "3", &err);   // x₀ = 3
    OcasHyperDual *y = ocas_dual_variable(shape, 1, "5", &err);   // x₁ = 5

    OcasHyperDual *x2 = ocas_dual_mul(x, x, &err);       // x², ∂/∂x = 6
    OcasHyperDual *xy = ocas_dual_mul(x, y, &err);       // xy,  ∂/∂x = 5
    OcasHyperDual *y2 = ocas_dual_mul(y, y, &err);       // y²,  ∂/∂y = 10
    OcasHyperDual *two_xy = ocas_dual_add(xy, xy, &err);
    OcasHyperDual *a = ocas_dual_add(x2, two_xy, &err);
    OcasHyperDual *f = ocas_dual_add(a, y2, &err);

    char *v = ocas_dual_value(f, &err);
    char *dfx = ocas_dual_deriv(f, 0, &err);
    char *dfy = ocas_dual_deriv(f, 1, &err);
    printf("f(3,5) = %s, ∂f/∂x = %s, ∂f/∂y = %s\n", v, dfx, dfy);

    ocas_string_free(v);
    ocas_string_free(dfx);
    ocas_string_free(dfy);
    ocas_hyperdual_free(f);
    ocas_hyperdual_free(a);
    ocas_hyperdual_free(two_xy);
    ocas_hyperdual_free(y2);
    ocas_hyperdual_free(xy);
    ocas_hyperdual_free(x2);
    ocas_hyperdual_free(y);
    ocas_hyperdual_free(x);
    ocas_dual_shape_free(shape);
    return 0;
}
```

> **注意**：双数域在 `ocas.hpp` 中**没有** C++ RAII 包装，必须直接使用上述 C API。

### C：数论

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // 素因数分解
    char *f = ocas_ntheory_factorint("360", &err);
    printf("360 = %s\n", f);  // 2:3,3:2,5:1
    ocas_string_free(f);

    // 素性测试
    printf("997 是素数？%s\n", ocas_ntheory_isprime("997", &err) ? "是" : "否");

    // 中国剩余定理
    char *r = ocas_ntheory_crt("3,5,7", "2,3,2", &err);
    printf("CRT: %s\n", r);  // 23,105
    ocas_string_free(r);

    // Euler totient
    char *t = ocas_ntheory_totient("100", &err);
    printf("φ(100) = %s\n", t);  // 40
    ocas_string_free(t);

    return 0;
}
```

---

## 参见

- [Rust API 参考](./rust.md) — Rust 公共 API 总览
- [Python API 参考](./python.md) — Python 绑定文档
- [架构](../architecture.md) — 项目架构概览
