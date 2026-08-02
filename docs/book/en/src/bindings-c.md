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

## ODE API

Since 0.20.1, `ocas-c` exposes ODE solving as string-in/string-out
functions. All returned strings are heap-allocated and must be released
with `ocas_string_free`.

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // Classify without solving.
    char *types = ocas_ode_classify(
        "Derivative(y(x), x) - y(x)", "y", "x", &err);
    printf("methods: %s\n", types);   /* LinearFirst,Separable,PowerSeries */

    // Solve (hint = NULL for auto classification).
    char *sol = ocas_ode_dsolve(
        "Derivative(y(x), x) - y(x)", "y", "x", NULL, &err);
    printf("solution: %s\n", sol);    /* y = C1*exp(x) */

    // Laplace IVP: y(0) = 1.
    char *ivp = ocas_ode_dsolve_ivp(
        "Derivative(y(x), x) - y(x)", "y", "x", "1", NULL, &err);
    printf("ivp: %s\n", ivp);          /* y = exp(x) */

    ocas_string_free(types);
    ocas_string_free(sol);
    ocas_string_free(ivp);
    return 0;
}
```

| Function | Purpose |
|---|---|
| `ocas_ode_classify` | Comma-separated applicable method names |
| `ocas_ode_dsolve` | Symbolic solution string (`"y = ..."` or `"unsolved"`) |
| `ocas_ode_dsolve_ivp` | Explicit IVP solution via Laplace transform |

## Numeric Integration API (Vegas)

Since 0.18.0, `ocas-c` provides Monte Carlo integration via opaque handles.

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // Create a 2-D Vegas integrator (n_samples=10000, iterations=10).
    ocas_OcasVegas* v = ocas_vegas_create(2, 10000, 10, &err);

    // Integrate using a callback.
    ocas_OcasIntegrateResult result;
    ocas_integrate_1d(my_fn, NULL, 0.0, 1.0, 10000, 10, &result);
    printf("integral = %f ± %f\n", result.integral, result.error);

    ocas_vegas_free(v);
    return 0;
}
```

| Function | Purpose |
|---|---|
| `ocas_vegas_create(n_dims, n_samples, iterations, &err)` | Create integrator |
| `ocas_vegas_integrate(vegas, fn, user_data, &err)` | Run integration |
| `ocas_vegas_result(vegas)` | Get `OcasIntegrateResult` |
| `ocas_integrate_1d(fn, user_data, a, b, n_samples, iterations, &result)` | 1-D convenience |

## Dual Numbers API (HyperDual)

Since 0.18.1, `ocas-c` provides forward-mode automatic differentiation
via opaque handles. Coefficients are strings (`"num"` or `"num/den"`).

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

| Function | Purpose |
|---|---|
| `ocas_dual_shape_new(n_vars, &err)` | First-order shape |
| `ocas_dual_variable(shape, i, coeff, &err)` | Variable with unit ε |
| `ocas_dual_constant(shape, coeff, &err)` | Constant (no ε) |
| `ocas_dual_value(hd, &err)` | Scalar value as string |
| `ocas_dual_deriv(hd, i, &err)` | i-th derivative as string |
| `ocas_dual_add/sub/mul/div(a, b, &err)` | Arithmetic |

## Tensor API

Since 0.18.1, `ocas-c` provides tensor creation, contraction, and
symmetrization. Index labels and positions are passed as arrays.

```c
#include <ocas.h>

int main(void) {
    int err = 0;

    // Create A^μ (rank 1, upper index "mu").
    const char* labels[] = {"mu"};
    int positions[] = {1};  /* 1 = upper */
    ocas_OcasTensor* A = ocas_tensor_create("A", labels, positions, 1, &err);

    int rank = ocas_tensor_rank(A, &err);  /* 1 */

    // Contract two tensors (returns a contraction result).
    ocas_OcasTensorContraction* c = ocas_tensor_contract(A, B, &err);
    /* Access scalar or free factors via c->scalar / c->product */

    ocas_tensor_contraction_free(c);
    ocas_tensor_free(A);
    return 0;
}
```

| Function | Purpose |
|---|---|
| `ocas_tensor_create(name, labels, positions, rank, &err)` | Create tensor |
| `ocas_tensor_rank(tensor, &err)` | Query rank |
| `ocas_tensor_symmetry(tensor, &err)` | Query symmetry (0=None, 1=Sym, 2=Anti) |
| `ocas_tensor_symmetrise_sign(tensor, &err)` | Antisymmetrization sign |
| `ocas_tensor_contract(a, b, &err)` | Contract matching indices |
| `ocas_tensor_to_string(tensor, &err)` | String representation |

## Algebraic Numbers API

Since 0.17.1, `ocas-c` provides algebraic number field operations via
opaque handles.

```c
#include <ocas.h>

int main(void) {
    int err = 0;

    // Q(√2): minimal polynomial x^2 - 2
    ocas_OcasAlgebraicField* field =
        ocas_algebraic_field_create("x^2 - 2", &err);
    int deg = ocas_algebraic_field_degree(field, &err);  /* 2 */

    // Create polynomial over Q(√2) and factor.
    ocas_OcasAlgebraicPoly* p =
        ocas_algebraic_poly_create(field, coeffs, n_coeffs, &err);
    ocas_OcasAlgebraicFactorArray factors;
    ocas_algebraic_poly_factor(p, &factors, &err);

    /* ... iterate factors ... */

    ocas_algebraic_factor_array_free(&factors);
    ocas_algebraic_poly_free(p);
    ocas_algebraic_field_free(field);
    return 0;
}
```

| Function | Purpose |
|---|---|
| `ocas_algebraic_field_create(min_poly, &err)` | Create from minimal polynomial |
| `ocas_algebraic_field_degree(field, &err)` | Extension degree |
| `ocas_algebraic_poly_create(field, coeffs, n, &err)` | Polynomial over the field |
| `ocas_algebraic_poly_factor(poly, &factors, &err)` | Factor (Trager's algorithm) |

## Gröbner Basis API

Since 0.23.0, `ocas-c` provides Gröbner basis computation and ideal operations
via coefficient arrays. Polynomials are represented as coefficient arrays with
explicit exponent matrices.

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // Polynomials: x² + y² - 1 and x - y in k[x,y]
    // Represented as coefficient arrays with exponent matrices.
    double coeffs1[] = {1.0, 1.0, -1.0};
    int exponents1[][2] = {{2, 0}, {0, 2}, {0, 0}};
    int n_terms1 = 3;
    int n_vars = 2;

    double coeffs2[] = {1.0, -1.0};
    int exponents2[][2] = {{1, 0}, {0, 1}};
    int n_terms2 = 2;

    // Compute Gröbner basis.
    OcasGroebnerBasis* gb = ocas_groebner_basis(
        coeffs1, exponents1, n_terms1,
        coeffs2, exponents2, n_terms2,
        n_vars, 0 /* algorithm: auto */, &err);

    int gb_len = ocas_groebner_basis_len(gb, &err);
    printf("GB has %d elements\n", gb_len);

    // Test membership.
    int contains = ocas_ideal_contains(
        coeffs1, exponents1, n_terms1,
        coeffs2, exponents2, n_terms2,
        n_vars, &err);
    printf("f in ideal: %d\n", contains);

    ocas_groebner_basis_free(gb);
    return 0;
}
```

| Function | Purpose |
|---|---|
| `ocas_groebner_basis(...)` | Compute Gröbner basis |
| `ocas_groebner_basis_len(gb, &err)` | Number of basis elements |
| `ocas_groebner_basis_free(gb)` | Release basis |
| `ocas_ideal_contains(...)` | Test ideal membership |
| `ocas_eliminate(...)` | Eliminate variables |
| `ocas_solve_system(...)` | Solve polynomial system |

Coefficients are `double` arrays; exponents are `int` matrices (row-major).
See [Gröbner Bases](./algorithms/groebner.md) for algorithmic details.
