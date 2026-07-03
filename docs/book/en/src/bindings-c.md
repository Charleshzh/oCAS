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

## Polynomial API

Since 0.11.1, `ocas-c` exposes bivariate polynomial objects as opaque
handles with factorization support over $\mathbb{Z}$ and $\mathbb{F}_p$.

### Integer polynomials (`OcasPolyZ`)

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err;
    // Create a bivariate integer polynomial from a string.
    OcasPolyZ* p = ocas_poly_z_create("x^2 + y + 1", &err);

    // Query total degree.
    printf("degree: %zu\n", ocas_poly_z_degree(p));

    // Factor into irreducible components.
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

### Finite-field polynomials (`OcasPolyFp`)

```c
// Create a polynomial over F_5 and factor it.
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

### Lifecycle

| Function | Purpose |
|---|---|
| `ocas_poly_z_create` / `ocas_poly_fp_create` | Create from string |
| `ocas_poly_z_clone` / `ocas_poly_fp_clone` | Deep copy |
| `ocas_poly_z_degree` / `ocas_poly_fp_degree` | Total degree |
| `ocas_poly_z_to_string` / `ocas_poly_fp_to_string` | Heap-allocated string (caller frees) |
| `ocas_poly_z_factor` / `ocas_poly_fp_factor` | Factor into irreducible components |
| `ocas_poly_z_free` / `ocas_poly_fp_free` | Release handle |
| `ocas_poly_factor_array_free` | Release factor array |

All polynomial functions are safe to call (no `unsafe` required). Passing
`NULL` to any function sets the error code and returns `NULL` / error.
