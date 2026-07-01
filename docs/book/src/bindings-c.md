# C/C++ API / C/C++ API

**English**

The `ocas-c` crate provides a stable C ABI (generated with `cbindgen`) for
expression lifecycle, calculus, and simplification, plus a C++ RAII wrapper
in `ocas-c/include/ocas.hpp`.

**中文**

`ocas-c` crate 提供稳定的 C ABI（由 `cbindgen` 生成），覆盖表达式生命周期、微积分与化简，并在 `ocas-c/include/ocas.hpp` 提供 C++ RAII 包装。

## Building / 构建

```bash
cargo build -p ocas-c --release
```

The shared library and the `ocas.h` / `ocas.hpp` headers are placed under
`ocas-c/include/`.

## C example / C 示例

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

## C++ RAII / C++ RAII

```cpp
#include <ocas.hpp>

int main() {
    ocas::Expression e("x^2 + 2*x + 1");
    auto d = e.diff("x");
    std::cout << d.to_string() << std::endl;   // 2*x + 2
    return 0;   // automatic cleanup
}
```

**English**

The C++ wrapper translates oCAS errors into `ocas::Error` exceptions and
manages arena-backed expressions via RAII, so manual `free` calls are
unnecessary.

**中文**

C++ 包装将 oCAS 错误转换为 `ocas::Error` 异常，并通过 RAII 管理 arena 后端表达式，无需手动调用 `free`。
