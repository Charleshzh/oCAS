# C/C++ API

The `ocas-c` crate provides a stable C ABI (generated with `cbindgen`) for
expression lifecycle, calculus, and simplification, plus a C++ RAII wrapper
in `ocas-c/include/ocas.hpp`.

## Building

```bash
cargo build -p ocas-c --release
```

The shared library and the `ocas.h` / `ocas.hpp` headers are placed under
`ocas-c/include/`.

## C example

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
    return 0;   // automatic cleanup
}
```

The C++ wrapper translates oCAS errors into `ocas::Error` exceptions and
manages arena-backed expressions via RAII, so manual `free` calls are
unnecessary.
