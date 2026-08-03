# C/C++ API Reference

> **⚠️ Important note: the `ocas.h` header is out of date**
>
> The cbindgen-generated `ocas-c/include/ocas.h` is missing the following 9 functions:
> - `ocas_expr_integrate_heuristic` (expression domain)
> - 8 Gröbner-domain exports (`ocas_groebner_basis`, `ocas_groebner_basis_free`, `ocas_groebner_basis_len`, `ocas_is_zero_dimensional`, `ocas_solve_polynomial_system`, `ocas_system_solution_count`, `ocas_system_solution_value`, `ocas_system_solution_free`)
>
> This document is based on the **actual exports** in `ocas-c/src/*.rs`, not the header file. When using the missing functions you must declare the prototypes yourself.

> **⚠️ Gröbner module caveats**
>
> The Gröbner module uses a **non-standard** `-1` error code (instead of the `OCAS_*` constants), and the `OcasGroebnerBasis` and `OcasSystemSolution` handles are **not marked** `#[repr(C)]` — they work correctly in the current Rust version, but ABI stability is not guaranteed by the compiler.

> **⚠️ C++ wrapper coverage**
>
> `ocas.hpp` covers only three domains (22 wrapper functions in total):
> - **Expressions**: the `ocas::Expression` class (8 members: 2 constructors + 6 methods)
> - **Number theory**: `ocas::ntheory::*` (11 inline functions)
> - **Tensors**: `ocas::tensor::*` (3 inline functions)
>
> The remaining 69 exports have no C++ RAII wrappers and must be called through the C API directly.

---

## General Conventions

### Error Handling

All fallible functions follow this pattern:

```c
int err = 0;
OcasExpr *e = ocas_expr_parse("x^2 + 1", &err);
if (e == NULL) {
    fprintf(stderr, "error %d: %s\n", err, ocas_error_last_message());
    ocas_error_clear();
    return 1;
}
```

- The `err_out` parameter may be `NULL` (error code ignored).
- Functions returning pointers: `NULL` on failure.
- Functions returning `int`: non-zero on failure (see each function's documentation for the exact meaning).
- Functions returning `size_t`: `0` on failure.
- Number theory functions returning `int` use the special sentinel values `-1` (`isprime`) or `-2` (`jacobi`/`mobius`/`liouville`) to signal errors.

### Memory Management

| Return type | Free function |
|---|---|
| `char*` (string) | `ocas_string_free()` |
| `OcasExpr*` | `ocas_expr_free()` |
| `OcasPolyZ*` | `ocas_poly_z_free()` |
| `OcasPolyFp*` | `ocas_poly_fp_free()` |
| `OcasPolyFactorArray` | `ocas_poly_factor_array_free()` (frees only the array structure, not the internal polynomial handles) |
| `OcasAlgebraicField*` | `ocas_algebraic_field_free()` |
| `OcasAlgebraicPoly*` | `ocas_algebraic_poly_free()` |
| `OcasAlgebraicFactorArray` | `ocas_algebraic_factor_array_free()` (as above, frees only the array) |
| `OcasVegas*` | `ocas_vegas_free()` |
| `OcasTensor*` | `ocas_tensor_free()` |
| `OcasTensorContraction` | `ocas_tensor_contraction_free()` (frees arrays/strings, not the internal tensor handles) |
| `OcasDualShape*` | `ocas_dual_shape_free()` |
| `OcasHyperDual*` | `ocas_hyperdual_free()` |
| `OcasGroebnerBasis*` | `ocas_groebner_basis_free()` |
| `OcasSystemSolution*` | `ocas_system_solution_free()` |

All `*_free` functions are safe on `NULL` arguments (no-op).

---

## 1. Errors/Utilities (3 functions)

### `ocas_version`

**Signature**:
```c
const char *ocas_version(void);
```

**Description**: Returns the oCAS version string.

**Return value**: A static string pointer valid for the lifetime of the program. **The caller must not free or modify it.**

**Example**:
```c
printf("oCAS version: %s\n", ocas_version());
// Output: oCAS version: 0.24.x
```

---

### `ocas_error_last_message`

**Signature**:
```c
const char *ocas_error_last_message(void);
```

**Description**: Returns the message of the last error on the calling thread, or `NULL` if there was no error.

**Return value**: A library-owned string pointer. **Must not be freed or modified.** It becomes invalid after calling any error-setting function or `ocas_error_clear()` on the same thread.

**Example**:
```c
ocas_expr_parse("invalid !!!", NULL);
const char *msg = ocas_error_last_message();
if (msg) printf("error: %s\n", msg);
ocas_error_clear();
```

---

### `ocas_error_clear`

**Signature**:
```c
void ocas_error_clear(void);
```

**Description**: Clears the last-error state on the calling thread.

**Example**:
```c
ocas_error_clear(); // reset the error state
```

---

## 2. Expressions (12 functions)

Expressions are the core type of oCAS. Each `OcasExpr*` handle owns a private arena, an `AtomArena`, and a normalized root `Atom`. Handles are movable; the heap address does not change over the handle's lifetime.

### `ocas_expr_parse`

**Signature**:
```c
OcasExpr *ocas_expr_parse(const char *input, int *err_out);
```

**Description**: Parses an input string into a new expression handle.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `input` | `const char*` | Null-terminated expression string, e.g. `"x^2 + 2*x + 1"` |
| `err_out` | `int*` | May be `NULL`; when non-NULL, the error code is written here |

**Return value**: `OcasExpr*` on success; `NULL` on failure (`OCAS_ERROR_PARSE`).

**Memory**: The caller is responsible for freeing with `ocas_expr_free()`.

**Example**:
```c
int err = 0;
OcasExpr *e = ocas_expr_parse("sin(x)^2 + cos(x)^2", &err);
assert(err == OCAS_OK && e != NULL);
ocas_expr_free(e);
```

---

### `ocas_expr_free`

**Signature**:
```c
void ocas_expr_free(OcasExpr *handle);
```

**Description**: Frees an expression handle. `NULL`-safe (no-op).

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `OcasExpr*` | Handle returned by `ocas_expr_parse` or other expression functions |

**Example**:
```c
ocas_expr_free(e);     // free
ocas_expr_free(NULL);  // no-op
```

---

### `ocas_expr_clone`

**Signature**:
```c
OcasExpr *ocas_expr_clone(const OcasExpr *handle, int *err_out);
```

**Description**: Clones the expression into a new arena. Returns an independent copy.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `const OcasExpr*` | Non-NULL expression handle |
| `err_out` | `int*` | May be `NULL` |

**Return value**: A new handle or `NULL` (`OCAS_ERROR_RUNTIME`).

**Example**:
```c
OcasExpr *copy = ocas_expr_clone(original, NULL);
// copy and original are independent; they do not affect each other
ocas_expr_free(copy);
```

---

### `ocas_expr_to_string`

**Signature**:
```c
char *ocas_expr_to_string(const OcasExpr *handle, int *err_out);
```

**Description**: Renders the expression as a null-terminated C string.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `const OcasExpr*` | Non-NULL expression handle |
| `err_out` | `int*` | May be `NULL` |

**Return value**: A heap-allocated string (free with `ocas_string_free()`), or `NULL`.

**Example**:
```c
char *s = ocas_expr_to_string(e, NULL);
printf("expression: %s\n", s);
ocas_string_free(s);
```

---

### `ocas_string_free`

**Signature**:
```c
void ocas_string_free(char *s);
```

**Description**: Frees a string allocated by `ocas_expr_to_string` or another function returning `char*`. `NULL`-safe.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `s` | `char*` | Heap-allocated string pointer |

---

### `ocas_expr_normalize`

**Signature**:
```c
int ocas_expr_normalize(OcasExpr *handle, int *err_out);
```

**Description**: Re-normalizes the expression handle in place.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `OcasExpr*` | Non-NULL expression handle (**modified**) |
| `err_out` | `int*` | May be `NULL` |

**Return value**: `OCAS_OK` (0) on success; non-zero on failure.

---

### `ocas_expr_diff`

**Signature**:
```c
OcasExpr *ocas_expr_diff(const OcasExpr *handle, const char *var, int *err_out);
```

**Description**: Differentiates the expression with respect to the variable `var`. Returns a new handle.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `const OcasExpr*` | Non-NULL expression handle |
| `var` | `const char*` | Variable name, e.g. `"x"` |
| `err_out` | `int*` | May be `NULL` |

**Return value**: The derivative expression handle or `NULL`.

**Example**:
```c
OcasExpr *e = ocas_expr_parse("x^3 + 2*x", NULL);
OcasExpr *de = ocas_expr_diff(e, "x", NULL);
char *s = ocas_expr_to_string(de, NULL);
printf("d/dx = %s\n", s);  // Output: 2 + (3*(x^2))
ocas_string_free(s);
ocas_expr_free(de);
ocas_expr_free(e);
```

---

### `ocas_expr_integrate`

**Signature**:
```c
OcasExpr *ocas_expr_integrate(const OcasExpr *handle, const char *var, int *err_out);
```

**Description**: Integrates the expression with respect to the variable `var`. If no symbolic solution is found, returns the unevaluated form `Integral(expr, var)`.

**Parameters**: Same as `ocas_expr_diff`.

**Return value**: The integral expression handle or `NULL`.

**Example**:
```c
OcasExpr *e = ocas_expr_parse("3*x^2", NULL);
OcasExpr *ie = ocas_expr_integrate(e, "x", NULL);
char *s = ocas_expr_to_string(ie, NULL);
printf("∫ = %s\n", s);  // Output: 3*(3^-1)*(x^3)
ocas_string_free(s);
ocas_expr_free(ie);
ocas_expr_free(e);
```

---

### `ocas_expr_integrate_heuristic`

**Signature**:
```c
OcasExpr *ocas_expr_integrate_heuristic(const OcasExpr *handle, const char *var, int *err_out);
```

**Description**: Integrates using heuristic techniques (integration by parts, trigonometric substitution, Weierstrass rationalization, Euler substitution). If no heuristic succeeds, returns the unevaluated form.

> **⚠️ Note**: This function is **not included in the `ocas.h` header**. Declare the prototype yourself:
> ```c
> extern OcasExpr *ocas_expr_integrate_heuristic(const OcasExpr *handle,
>                                                 const char *var,
>                                                 int *err_out);
> ```

**Parameters**: Same as `ocas_expr_diff`.

**Return value**: The integral expression handle or `NULL`.

---

### `ocas_expr_taylor`

**Signature**:
```c
OcasExpr *ocas_expr_taylor(const OcasExpr *handle, const char *var,
                           const OcasExpr *point, uint32_t order, int *err_out);
```

**Description**: Computes the Taylor expansion of the expression in `var` of order `order` at the point `point`.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `const OcasExpr*` | Non-NULL expression handle |
| `var` | `const char*` | Expansion variable |
| `point` | `const OcasExpr*` | Expansion point (expression handle) |
| `order` | `uint32_t` | Expansion order |
| `err_out` | `int*` | May be `NULL` |

**Return value**: The Taylor expansion handle or `NULL`.

**Example**:
```c
OcasExpr *f = ocas_expr_parse("exp(x)", NULL);
OcasExpr *zero = ocas_expr_parse("0", NULL);
OcasExpr *t = ocas_expr_taylor(f, "x", zero, 4, NULL);
char *s = ocas_expr_to_string(t, NULL);
printf("Taylor = %s\n", s);  // Output: 1 + x + ((2^-1)*(x^2)) + ((6^-1)*(x^3)) + ((24^-1)*(x^4))
ocas_string_free(s);
ocas_expr_free(t);
ocas_expr_free(zero);
ocas_expr_free(f);
```

---

### `ocas_expr_simplify`

**Signature**:
```c
OcasExpr *ocas_expr_simplify(const OcasExpr *handle, int *err_out);
```

**Description**: Simplifies the expression using the default rule set. Returns a new handle.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `const OcasExpr*` | Non-NULL expression handle |
| `err_out` | `int*` | May be `NULL` |

**Return value**: The simplified expression handle or `NULL`.

**Example**:
```c
OcasExpr *e = ocas_expr_parse("x + x", NULL);
OcasExpr *simplified = ocas_expr_simplify(e, NULL);
char *s = ocas_expr_to_string(simplified, NULL);
printf("simplified = %s\n", s);  // Output: 2*x
ocas_string_free(s);
ocas_expr_free(simplified);
ocas_expr_free(e);
```

---

### `ocas_expr_substitute`

**Signature**:
```c
OcasExpr *ocas_expr_substitute(const OcasExpr *handle, const char *var,
                                const OcasExpr *replacement, int *err_out);
```

**Description**: Replaces every occurrence of `var` in the expression with `replacement`.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `handle` | `const OcasExpr*` | Non-NULL expression handle |
| `var` | `const char*` | Name of the variable to replace |
| `replacement` | `const OcasExpr*` | Replacement expression (non-NULL) |
| `err_out` | `int*` | May be `NULL` |

**Return value**: The substituted expression handle or `NULL`.

**Example**:
```c
OcasExpr *e = ocas_expr_parse("x^2 + y", NULL);
OcasExpr *val = ocas_expr_parse("3", NULL);
OcasExpr *result = ocas_expr_substitute(e, "x", val, NULL);
char *s = ocas_expr_to_string(result, NULL);
printf("substituted = %s\n", s);  // Output: y + (3^2)
ocas_string_free(s);
ocas_expr_free(result);
ocas_expr_free(val);
ocas_expr_free(e);
```

---

### C++ wrapper: `ocas::Expression`

`ocas.hpp` provides a RAII wrapper class that manages handle lifetimes automatically:

```cpp
#include <ocas.h>
#include <ocas.hpp>

try {
    ocas::Expression e("x^2 + 1");
    ocas::Expression d = e.diff("x");
    std::cout << "d/dx = " << d.to_string() << std::endl;
    // all handles freed automatically
} catch (const ocas::Error& ex) {
    std::cerr << "error: " << ex.what() << std::endl;
}
```

| C++ method | Corresponding C function |
|---|---|
| `Expression(input)` | `ocas_expr_parse` |
| `Expression(other)` | `ocas_expr_clone` |
| `to_string()` | `ocas_expr_to_string` + `ocas_string_free` |
| `diff(var)` | `ocas_expr_diff` |
| `integrate(var)` | `ocas_expr_integrate` |
| `simplify()` | `ocas_expr_simplify` |
| `substitute(var, rep)` | `ocas_expr_substitute` |
| `raw()` | Get the raw `OcasExpr*` (does not transfer ownership) |

---

## 3. Integer Polynomials (6 functions)

Handles bivariate (`x`, `y`) integer polynomials ℤ[x,y]. Parsed from strings; factorization supported.

### `ocas_poly_z_create`

**Signature**:
```c
OcasPolyZ *ocas_poly_z_create(const char *input, int *err);
```

**Description**: Creates a bivariate integer polynomial from a string. The input may contain the variables `x`, `y`, integer coefficients, addition, multiplication, and non-negative integer powers.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `input` | `const char*` | Polynomial string, e.g. `"x^2*y + 3*x - 1"` |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasPolyZ*` handle or `NULL` (`OCAS_ERROR_PARSE`).

**Memory**: The caller frees with `ocas_poly_z_free()`.

**Example**:
```c
int err = 0;
OcasPolyZ *p = ocas_poly_z_create("x^2 + y^2 - 1", &err);
assert(p != NULL);
ocas_poly_z_free(p);
```

---

### `ocas_poly_z_free`

**Signature**:
```c
void ocas_poly_z_free(OcasPolyZ *poly);
```

**Description**: Frees an integer polynomial handle. `NULL`-safe.

---

### `ocas_poly_z_clone`

**Signature**:
```c
OcasPolyZ *ocas_poly_z_clone(const OcasPolyZ *poly);
```

**Description**: Clones an integer polynomial. Returns a new handle (caller frees it).

**Return value**: A new `OcasPolyZ*` or `NULL` (NULL input or out of memory).

---

### `ocas_poly_z_degree`

**Signature**:
```c
size_t ocas_poly_z_degree(const OcasPolyZ *poly);
```

**Description**: Returns the **total degree** of the polynomial. The zero polynomial returns `0`. A `NULL` handle returns `0`.

---

### `ocas_poly_z_to_string`

**Signature**:
```c
char *ocas_poly_z_to_string(const OcasPolyZ *poly, int *err);
```

**Description**: Returns a heap-allocated string representation of the polynomial. Free with `ocas_string_free()`.

**Return value**: The string or `NULL`.

**Example**:
```c
char *s = ocas_poly_z_to_string(p, NULL);
printf("polynomial: %s\n", s);
ocas_string_free(s);
```

---

### `ocas_poly_z_factor`

**Signature**:
```c
int ocas_poly_z_factor(const OcasPolyZ *poly, OcasPolyFactorArray *out, int *err);
```

**Description**: Factors a bivariate integer polynomial.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `poly` | `const OcasPolyZ*` | Non-NULL polynomial handle |
| `out` | `OcasPolyFactorArray*` | Output struct (filled on success) |
| `err` | `int*` | May be `NULL` |

**Return value**: `OCAS_OK` (0) on success; non-zero on failure (`out` unchanged).

**Memory**: On success, `out.factors[i].poly` is a `void*` that must be cast to `OcasPolyZ*` and freed with `ocas_poly_z_free()`. The array itself is freed with `ocas_poly_factor_array_free()`.

**Example**:
```c
OcasPolyZ *p = ocas_poly_z_create("x^2 - y^2", NULL);
OcasPolyFactorArray arr = {0};
if (ocas_poly_z_factor(p, &arr, NULL) == OCAS_OK) {
    for (size_t i = 0; i < arr.len; i++) {
        OcasPolyZ *f = (OcasPolyZ *)arr.factors[i].poly;
        char *s = ocas_poly_z_to_string(f, NULL);
        printf("factor: %s (multiplicity %zu)\n", s, arr.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_z_free(f);
    }
    ocas_poly_factor_array_free(&arr);
}
ocas_poly_z_free(p);
```

---

## 4. Finite Field Polynomials (6 functions)

Handles bivariate polynomials $𝔽_p[x,y]$.

### `ocas_poly_fp_create`

**Signature**:
```c
OcasPolyFp *ocas_poly_fp_create(const char *input, const char *prime, int *err);
```

**Description**: Creates a bivariate polynomial over the prime field $𝔽_p$. Coefficients in the string are automatically reduced modulo $p$.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `input` | `const char*` | Polynomial string |
| `prime` | `const char*` | Decimal string of the prime $p$ |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasPolyFp*` handle or `NULL`.

**Memory**: The caller frees with `ocas_poly_fp_free()`.

**Example**:
```c
OcasPolyFp *p = ocas_poly_fp_create("x^2 + 2*x + 1", "7", NULL);
// coefficients in 𝔽₇: x² + 2x + 1
ocas_poly_fp_free(p);
```

---

### `ocas_poly_fp_free`

**Signature**:
```c
void ocas_poly_fp_free(OcasPolyFp *poly);
```

**Description**: Frees a finite field polynomial handle. `NULL`-safe.

---

### `ocas_poly_fp_clone`

**Signature**:
```c
OcasPolyFp *ocas_poly_fp_clone(const OcasPolyFp *poly);
```

**Description**: Clones a finite field polynomial. Returns a new handle.

---

### `ocas_poly_fp_degree`

**Signature**:
```c
size_t ocas_poly_fp_degree(const OcasPolyFp *poly);
```

**Description**: Returns the total degree of the polynomial. `NULL` returns `0`.

---

### `ocas_poly_fp_to_string`

**Signature**:
```c
char *ocas_poly_fp_to_string(const OcasPolyFp *poly, int *err);
```

**Description**: Returns a heap-allocated string of the polynomial. Free with `ocas_string_free()`.

---

### `ocas_poly_fp_factor`

**Signature**:
```c
int ocas_poly_fp_factor(const OcasPolyFp *poly, OcasPolyFactorArray *out, int *err);
```

**Description**: Factors a finite field polynomial. Return value and memory management are the same as `ocas_poly_z_factor`, except that `out.factors[i].poly` must be cast to `OcasPolyFp*`.

**Example**:
```c
OcasPolyFp *p = ocas_poly_fp_create("x^2 - 1", "5", NULL);
OcasPolyFactorArray arr = {0};
if (ocas_poly_fp_factor(p, &arr, NULL) == OCAS_OK) {
    for (size_t i = 0; i < arr.len; i++) {
        OcasPolyFp *f = (OcasPolyFp *)arr.factors[i].poly;
        char *s = ocas_poly_fp_to_string(f, NULL);
        printf("factor: %s (multiplicity %zu)\n", s, arr.factors[i].multiplicity);
        ocas_string_free(s);
        ocas_poly_fp_free(f);
    }
    ocas_poly_factor_array_free(&arr);
}
ocas_poly_fp_free(p);
```

---

## 5. Factor Arrays (1 function)

### `ocas_poly_factor_array_free`

**Signature**:
```c
void ocas_poly_factor_array_free(OcasPolyFactorArray *arr);
```

**Description**: Frees the factor array structure itself and the internal `OcasPolyFactor` objects, but **does not free** the polynomial handles of the individual factors.

> **⚠️ Memory pitfall**: you must free each `arr->factors[i].poly` with `ocas_poly_z_free()` / `ocas_poly_fp_free()` *before* calling this function. Doing it in the wrong order causes a use-after-free.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `arr` | `OcasPolyFactorArray*` | Array filled by `ocas_poly_z_factor` or `ocas_poly_fp_factor` |

---

## 6. Algebraic Number Fields (9 functions)

Handles polynomials over algebraic number fields $ℚ(α)$. Supports the Trager factorization algorithm.

### Data Formats

**Minimal polynomial**: a univariate string; the variable must be `x`, e.g. `"x^2 - 2"` denotes $ℚ(\sqrt{2})$.

**Coefficient lists**: semicolon-separated (constant term first); each coefficient is a comma-separated list of rationals (an $α$-polynomial, ascending order). For example, over $ℚ(\sqrt{2})$:
- `"-2;0;1"` — $x^2 - 2$ (all coefficients in the base field)
- `"0;0;0,1"` — $x^2 - α$ (the $x^2$ coefficient is $0 + 1 \cdot α$)

A single rational may omit the comma.

### `ocas_algebraic_field_create`

**Signature**:
```c
OcasAlgebraicField *ocas_algebraic_field_create(const char *min_poly, int *err);
```

**Description**: Creates an algebraic number field from a monic minimal polynomial.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `min_poly` | `const char*` | Minimal polynomial string, e.g. `"x^2 - 2"` |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasAlgebraicField*` handle or `NULL`.

**Memory**: The caller frees with `ocas_algebraic_field_free()`.

**Example**:
```c
OcasAlgebraicField *q_sqrt2 = ocas_algebraic_field_create("x^2 - 2", NULL);
printf("extension degree: %zu\n", ocas_algebraic_field_degree(q_sqrt2));
// Output: extension degree: 2
ocas_algebraic_field_free(q_sqrt2);
```

---

### `ocas_algebraic_field_free`

**Signature**:
```c
void ocas_algebraic_field_free(OcasAlgebraicField *field);
```

**Description**: Frees an algebraic number field handle. `NULL`-safe.

---

### `ocas_algebraic_field_degree`

**Signature**:
```c
size_t ocas_algebraic_field_degree(const OcasAlgebraicField *field);
```

**Description**: Returns the extension degree $\deg(m) = [ℚ(α):ℚ]$. `NULL` returns `0`.

---

### `ocas_algebraic_poly_create`

**Signature**:
```c
OcasAlgebraicPoly *ocas_algebraic_poly_create(const OcasAlgebraicField *field,
                                              const char *coeffs, int *err);
```

**Description**: Creates a polynomial over the algebraic number field.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `field` | `const OcasAlgebraicField*` | Non-NULL field handle |
| `coeffs` | `const char*` | Coefficient list string (see format above) |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasAlgebraicPoly*` handle or `NULL`.

**Example**:
```c
OcasAlgebraicField *fld = ocas_algebraic_field_create("x^2 - 2", NULL);
// Create the polynomial x² - α (i.e. x² - √2)
OcasAlgebraicPoly *p = ocas_algebraic_poly_create(fld, "0;0;0,1", NULL);
printf("degree: %zu\n", ocas_algebraic_poly_degree(p));
ocas_algebraic_poly_free(p);
ocas_algebraic_field_free(fld);
```

---

### `ocas_algebraic_poly_free`

**Signature**:
```c
void ocas_algebraic_poly_free(OcasAlgebraicPoly *poly);
```

**Description**: Frees an algebraic field polynomial handle. `NULL`-safe.

---

### `ocas_algebraic_poly_degree`

**Signature**:
```c
size_t ocas_algebraic_poly_degree(const OcasAlgebraicPoly *poly);
```

**Description**: Returns the degree of the polynomial. The zero polynomial returns `0`.

---

### `ocas_algebraic_poly_to_string`

**Signature**:
```c
char *ocas_algebraic_poly_to_string(const OcasAlgebraicPoly *poly, int *err);
```

**Description**: Returns a heap-allocated string of the polynomial. The format is `[c0] + [c1]*x + [c2]*x^2 + ...`, where each `[ci]` is a comma-separated list of $α$-polynomial rationals. Free with `ocas_string_free()`.

---

### `ocas_algebraic_poly_factor`

**Signature**:
```c
int ocas_algebraic_poly_factor(const OcasAlgebraicPoly *poly,
                               OcasAlgebraicFactorArray *out, int *err);
```

**Description**: Factors a polynomial over an algebraic number field using the Trager algorithm.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `poly` | `const OcasAlgebraicPoly*` | Non-NULL polynomial handle |
| `out` | `OcasAlgebraicFactorArray*` | Output array (filled on success) |
| `err` | `int*` | May be `NULL` |

**Return value**: `OCAS_OK` (0) on success.

**Memory**: `out->factors[i].poly` (a `void*`) must be cast to `OcasAlgebraicPoly*` and freed with `ocas_algebraic_poly_free()`. The array is freed with `ocas_algebraic_factor_array_free()`.

**Example**:
```c
OcasAlgebraicField *fld = ocas_algebraic_field_create("x^2 - 2", NULL);
OcasAlgebraicPoly *p = ocas_algebraic_poly_create(fld, "-2;0;1", NULL);
OcasAlgebraicFactorArray arr = {0};
if (ocas_algebraic_poly_factor(p, &arr, NULL) == OCAS_OK) {
    for (size_t i = 0; i < arr.len; i++) {
        OcasAlgebraicPoly *f = (OcasAlgebraicPoly *)arr.factors[i].poly;
        char *s = ocas_algebraic_poly_to_string(f, NULL);
        printf("factor: %s (multiplicity %zu)\n", s, arr.factors[i].multiplicity);
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

**Signature**:
```c
void ocas_algebraic_factor_array_free(OcasAlgebraicFactorArray *arr);
```

**Description**: Frees the storage of the factor array. **Does not free** the polynomial handles of the individual factors — free those separately first.

---

## 7. Numeric Integration (6 functions)

Vegas adaptive Monte Carlo integrator.

### Integrand function signature

```c
typedef double (*ocas_integrand_t)(double x, void *user_data);
```

`user_data` is passed through from the caller to the integrand unchanged.

> **⚠️ Note**: the current C API passes only the **first coordinate** of the sampled point to the integrand (a 1-D-friendly interface). Even if the integrator is created with `n_dims > 1`, the `x` argument of the integrand represents the first dimension only; a full multi-dimensional `const double*`/`size_t` signature is planned for a future release.

### `ocas_vegas_create`

**Signature**:
```c
OcasVegas *ocas_vegas_create(size_t n_dims, const OcasVegasOptions *opts, int *err);
```

**Description**: Creates an `n_dims`-dimensional Vegas integrator.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `n_dims` | `size_t` | Number of integration dimensions |
| `opts` | `const OcasVegasOptions*` | May be `NULL` (uses defaults) |
| `err` | `int*` | May be `NULL` |

**`OcasVegasOptions` fields**:
| Field | Type | Default | Description |
|---|---|---|---|
| `n_bins` | `size_t` | 64 | Number of bins per dimension |
| `n_samples` | `size_t` | 10000 | Samples per iteration |
| `iterations` | `size_t` | 10 | Number of adaptive iterations |
| `learning_rate` | `double` | 1.5 | Grid smoothing/learning rate |
| `seed` | `uint64_t` | 0x0C45 | Random number seed |

**Return value**: `OcasVegas*` handle or `NULL`.

**Memory**: The caller frees with `ocas_vegas_free()`.

**Example**:
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

**Signature**:
```c
void ocas_vegas_free(OcasVegas *v);
```

**Description**: Frees a Vegas integrator handle. `NULL`-safe.

---

### `ocas_vegas_integrate`

**Signature**:
```c
OcasIntegrateResult ocas_vegas_integrate(OcasVegas *v, ocas_integrand_t f,
                                         void *user_data, int *err);
```

**Description**: Integrates `f` over the unit hypercube (the integrand receives only the first coordinate; see the note above).

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `v` | `OcasVegas*` | Non-NULL integrator handle |
| `f` | `ocas_integrand_t` | Integrand function |
| `user_data` | `void*` | User data passed to the integrand |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasIntegrateResult { double integral; double error; }`.

**Example**:
```c
static double f(double x, void *ud) { return x * x; }

OcasIntegrateResult r = ocas_vegas_integrate(v, f, NULL, &err);
printf("∫₀¹ x² dx ≈ %g ± %g\n", r.integral, r.error);
// Output (sample values; they depend on the seed and sample count): ∫₀¹ x² dx ≈ 0.333 ± 0.001
```

---

### `ocas_vegas_result`

**Signature**:
```c
OcasIntegrateResult ocas_vegas_result(const OcasVegas *v);
```

**Description**: Returns the cumulative estimate and error of the most recent `ocas_vegas_integrate` call.

---

### `ocas_vegas_iterations`

**Signature**:
```c
size_t ocas_vegas_iterations(const OcasVegas *v);
```

**Description**: Returns the number of completed iterations. `NULL` returns `0`.

---

### `ocas_integrate_1d`

**Signature**:
```c
OcasIntegrateResult ocas_integrate_1d(ocas_integrand_t f, void *user_data,
                                      double a, double b,
                                      const OcasVegasOptions *opts, int *err);
```

**Description**: Convenience function for one-dimensional numeric integration — integrates over $[a, b]$ using Vegas.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `f` | `ocas_integrand_t` | Integrand |
| `user_data` | `void*` | User data |
| `a` | `double` | Lower integration limit |
| `b` | `double` | Upper integration limit |
| `opts` | `const OcasVegasOptions*` | May be `NULL` |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasIntegrateResult`.

**Example**:
```c
OcasIntegrateResult r = ocas_integrate_1d(f, NULL, 0.0, 1.0, NULL, &err);
printf("result: %g ± %g\n", r.integral, r.error);
```

---

## 8. ODE (3 functions)

Symbolic ordinary differential equation solving. All inputs and outputs are strings.

### `ocas_ode_classify`

**Signature**:
```c
char *ocas_ode_classify(const char *equation, const char *func,
                        const char *var, int *err);
```

**Description**: Classifies an ODE and returns a comma-separated list of applicable solution method names.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `equation` | `const char*` | Expression equal to zero, e.g. `"Derivative(y(x), x) - y(x)"` |
| `func` | `const char*` | Unknown function name, e.g. `"y"` |
| `var` | `const char*` | Independent variable, e.g. `"x"` |
| `err` | `int*` | May be `NULL` |

**Return value**: A heap-allocated string such as `"LinearFirst,PowerSeries"`; free with `ocas_string_free()`. `NULL` on failure.

**Example**:
```c
char *methods = ocas_ode_classify("Derivative(y(x),x) - y(x)", "y", "x", NULL);
printf("available methods: %s\n", methods);
// Output: available methods: LinearFirst,PowerSeries
ocas_string_free(methods);
```

---

### `ocas_ode_dsolve`

**Signature**:
```c
char *ocas_ode_dsolve(const char *equation, const char *func,
                      const char *var, const char *hint, int *err);
```

**Description**: Solves an ODE symbolically.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `equation` | `const char*` | ODE expression |
| `func` | `const char*` | Unknown function name |
| `var` | `const char*` | Independent variable |
| `hint` | `const char*` | May be `NULL` (automatic classification); or one of the method names returned by `ocas_ode_classify` |
| `err` | `int*` | May be `NULL` |

**Return value**: A solution string such as `"y = C1*exp(x)"` or `"unsolved"`. Free with `ocas_string_free()`.

**Example**:
```c
char *sol = ocas_ode_dsolve("Derivative(y(x),x) - y(x)", "y", "x", NULL, NULL);
printf("solution: %s\n", sol);  // Output: solution: y = C1*exp(x)
ocas_string_free(sol);
```

---

### `ocas_ode_dsolve_ivp`

**Signature**:
```c
char *ocas_ode_dsolve_ivp(const char *equation, const char *func,
                          const char *var, const char *y0,
                          const char *y1, int *err);
```

**Description**: Solves an initial value problem for a first- or second-order linear constant-coefficient ODE via the Laplace transform.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `equation` | `const char*` | ODE expression |
| `func` | `const char*` | Unknown function name |
| `var` | `const char*` | Independent variable |
| `y0` | `const char*` | String expression for $y(0)$, e.g. `"1"` |
| `y1` | `const char*` | $y'(0)$ (may be `NULL`; only needed for second order) |
| `err` | `int*` | May be `NULL` |

**Return value**: An explicit solution string without free constants. Free with `ocas_string_free()`.

**Example**:
```c
// y'' + y = 0, y(0) = 1, y'(0) = 0 → y = cos(x)
char *sol = ocas_ode_dsolve_ivp(
    "Derivative(y(x),x,2) + y(x)", "y", "x", "1", "0", NULL);
printf("IVP solution: %s\n", sol);
ocas_string_free(sol);
```

---

## 9. Number Theory (11 functions)

Arbitrary-precision integers are passed as decimal strings. String results are freed by the caller with `ocas_string_free()`.

### `ocas_ntheory_factorint`

**Signature**:
```c
char *ocas_ntheory_factorint(const char *n, int *err_out);
```

**Description**: Computes the prime factorization of $|n|$.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `n` | `const char*` | Decimal integer string |
| `err_out` | `int*` | May be `NULL` |

**Return value**: `"p1:e1,p2:e2,..."` (ascending); negative numbers start with `"-1:1"`. `NULL` on failure.

**Example**:
```c
char *f = ocas_ntheory_factorint("360", NULL);
printf("360 = %s\n", f);  // Output: 360 = 2:3,3:2,5:1
ocas_string_free(f);

f = ocas_ntheory_factorint("-12", NULL);
printf("-12 = %s\n", f);  // Output: -12 = -1:1,2:2,3:1
ocas_string_free(f);
```

---

### `ocas_ntheory_isprime`

**Signature**:
```c
int ocas_ntheory_isprime(const char *n, int *err_out);
```

**Description**: BPSW probable primality test.

**Return value**: `1` = (very likely) prime, `0` = composite, `-1` = error.

**Example**:
```c
int r = ocas_ntheory_isprime("997", NULL);
printf("is 997 prime? %s\n", r == 1 ? "yes" : "no");
// Output: is 997 prime? yes
```

---

### `ocas_ntheory_nextprime`

**Signature**:
```c
char *ocas_ntheory_nextprime(const char *n, int *err_out);
```

**Description**: Returns the smallest prime strictly greater than `n`.

**Return value**: A decimal string or `NULL`.

**Example**:
```c
char *p = ocas_ntheory_nextprime("100", NULL);
printf("next prime after 100: %s\n", p);  // Output: 101
ocas_string_free(p);
```

---

### `ocas_ntheory_discrete_log`

**Signature**:
```c
char *ocas_ntheory_discrete_log(const char *p, const char *base,
                                const char *target, int *err_out);
```

**Description**: Solves $\text{base}^x \equiv \text{target} \pmod{p}$. Uses Pohlig–Hellman for prime $p$, BSGS otherwise.

**Return value**: A decimal string of the logarithm, or `NULL` (sets `OCAS_ERROR_RUNTIME` if no solution exists).

**Example**:
```c
char *x = ocas_ntheory_discrete_log("7", "3", "5", NULL);
// 3^x ≡ 5 (mod 7) → x = 5
printf("discrete log: %s\n", x);
ocas_string_free(x);
```

---

### `ocas_ntheory_crt`

**Signature**:
```c
char *ocas_ntheory_crt(const char *moduli, const char *residues, int *err_out);
```

**Description**: Chinese remainder theorem. The moduli need not be coprime; returns `NULL` if inconsistent.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `moduli` | `const char*` | Comma-separated moduli, e.g. `"3,5,7"` |
| `residues` | `const char*` | Comma-separated residues, e.g. `"2,3,2"` |

**Return value**: `"r,m"` (with $r \equiv \text{residues}[i] \pmod{\text{moduli}[i]}$), or `NULL`.

**Example**:
```c
char *r = ocas_ntheory_crt("3,5,7", "2,3,2", NULL);
printf("CRT: %s\n", r);  // Output: CRT: 23,105
// 23 ≡ 2 (mod 3), 23 ≡ 3 (mod 5), 23 ≡ 2 (mod 7)
ocas_string_free(r);
```

---

### `ocas_ntheory_jacobi`

**Signature**:
```c
int ocas_ntheory_jacobi(const char *a, const char *n, int *err_out);
```

**Description**: Computes the Jacobi symbol $(a/n)$; `n` must be a positive odd integer.

**Return value**: `-1`, `0`, or `1`; `-2` for invalid input.

---

### `ocas_ntheory_totient`

**Signature**:
```c
char *ocas_ntheory_totient(const char *n, int *err_out);
```

**Description**: Euler totient function $\varphi(n)$.

**Return value**: A decimal string or `NULL`.

**Example**:
```c
char *t = ocas_ntheory_totient("12", NULL);
printf("φ(12) = %s\n", t);  // Output: φ(12) = 4
ocas_string_free(t);
```

---

### `ocas_ntheory_mobius`

**Signature**:
```c
int ocas_ntheory_mobius(const char *n, int *err_out);
```

**Description**: Möbius function $\mu(n)$.

**Return value**: `-1`, `0`, or `1`; `-2` for errors.

---

### `ocas_ntheory_divisor_count`

**Signature**:
```c
char *ocas_ntheory_divisor_count(const char *n, int *err_out);
```

**Description**: Number of positive divisors $\tau(n)$.

**Return value**: A decimal string or `NULL`.

**Example**:
```c
char *d = ocas_ntheory_divisor_count("12", NULL);
printf("τ(12) = %s\n", d);  // Output: τ(12) = 6
// divisors of 12: 1,2,3,4,6,12
ocas_string_free(d);
```

---

### `ocas_ntheory_divisor_sigma`

**Signature**:
```c
char *ocas_ntheory_divisor_sigma(const char *n, uint32_t k, int *err_out);
```

**Description**: Sum of the $k$-th powers of the positive divisors, $\sigma_k(n)$.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `n` | `const char*` | Positive integer |
| `k` | `uint32_t` | Power |
| `err_out` | `int*` | May be `NULL` |

**Return value**: A decimal string or `NULL`.

**Example**:
```c
char *s = ocas_ntheory_divisor_sigma("6", 1, NULL);
printf("σ₁(6) = %s\n", s);  // Output: σ₁(6) = 12
// 1+2+3+6 = 12
ocas_string_free(s);
```

---

### `ocas_ntheory_liouville`

**Signature**:
```c
int ocas_ntheory_liouville(const char *n, int *err_out);
```

**Description**: Liouville function $\lambda(n) = (-1)^{\Omega(n)}$.

**Return value**: `-1`, `0`, or `1`; `-2` for errors.

---

### C++ wrapper: `ocas::ntheory::*`

All number theory functions have C++ inline wrappers that manage string memory automatically:

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

| C++ function | C function | Description |
|---|---|---|
| `factorint(n)` | `ocas_ntheory_factorint` | Prime factorization |
| `isprime(n)` | `ocas_ntheory_isprime` | Returns `bool`; throws on error |
| `nextprime(n)` | `ocas_ntheory_nextprime` | Next prime |
| `discrete_log(p, base, target)` | `ocas_ntheory_discrete_log` | Discrete logarithm; throws if no solution |
| `crt(moduli, residues)` | `ocas_ntheory_crt` | Chinese remainder theorem |
| `jacobi(a, n)` | `ocas_ntheory_jacobi` | Jacobi symbol; throws on invalid input |
| `totient(n)` | `ocas_ntheory_totient` | Euler totient |
| `mobius(n)` | `ocas_ntheory_mobius` | Möbius function; throws on invalid input |
| `divisor_count(n)` | `ocas_ntheory_divisor_count` | Number of divisors |
| `divisor_sigma(n, k=1)` | `ocas_ntheory_divisor_sigma` | Sum of divisor powers |
| `liouville(n)` | `ocas_ntheory_liouville` | Liouville function; throws on invalid input |

All wrappers throw an `ocas::Error` exception when the underlying C API reports an error.

---

## 10. Tensors (12 functions)

Named-index tensors supporting index contraction and symmetry.

### Slot string format

Slots are semicolon-separated; each entry is `label,position`, where `position` is `"upper"` or `"lower"` (aliases: `up`/`down`/`contravariant`/`covariant`).

For example: `"i,upper;j,lower"` denotes two slots $i^j{}_{}$.

### `ocas_tensor_create`

**Signature**:
```c
OcasTensor *ocas_tensor_create(const char *name, const char *slots,
                               const char *symmetry, int *err);
```

**Description**: Creates a tensor.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `name` | `const char*` | Tensor name, e.g. `"T"` |
| `slots` | `const char*` | Slot string |
| `symmetry` | `const char*` | May be `NULL` (= `"none"`); or `"none"`/`"symmetric"`/`"antisymmetric"` |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasTensor*` handle or `NULL`.

**Example**:
```c
OcasTensor *T = ocas_tensor_create("T", "i,upper;j,lower", "symmetric", NULL);
printf("rank: %zu\n", ocas_tensor_rank(T));  // Output: rank: 2
ocas_tensor_free(T);
```

---

### `ocas_tensor_free`

**Signature**:
```c
void ocas_tensor_free(OcasTensor *t);
```

**Description**: Frees a tensor handle. `NULL`-safe.

---

### `ocas_tensor_name`

**Signature**:
```c
char *ocas_tensor_name(const OcasTensor *t, int *err);
```

**Description**: Returns a heap-allocated string of the tensor name. Free with `ocas_string_free()`.

---

### `ocas_tensor_rank`

**Signature**:
```c
size_t ocas_tensor_rank(const OcasTensor *t);
```

**Description**: Returns the rank of the tensor (number of slots). `NULL` returns `0`.

---

### `ocas_tensor_symmetry`

**Signature**:
```c
int ocas_tensor_symmetry(const OcasTensor *t);
```

**Description**: Returns the symmetry code.

**Return value**: `0` = none, `1` = symmetric, `2` = antisymmetric, `-1` = null handle.

---

### `ocas_tensor_to_string`

**Signature**:
```c
char *ocas_tensor_to_string(const OcasTensor *t, int *err);
```

**Description**: Returns the string representation `name(slot, slot, ...)` of the tensor. Free with `ocas_string_free()`.

---

### `ocas_tensor_symmetrise_sign`

**Signature**:
```c
int64_t ocas_tensor_symmetrise_sign(const OcasTensor *t);
```

**Description**: Returns the antisymmetrization sign ($+1$ or $-1$). `NULL` returns `0`.

---

### `ocas_tensor_contract`

**Signature**:
```c
int ocas_tensor_contract(const OcasTensor *a, const OcasTensor *b,
                         OcasTensorContraction *out, int *err);
```

**Description**: Contracts over the shared dummy indices (identical labels with opposite variance) of two tensors.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `a` | `const OcasTensor*` | Non-NULL tensor |
| `b` | `const OcasTensor*` | Non-NULL tensor |
| `out` | `OcasTensorContraction*` | Output struct |
| `err` | `int*` | May be `NULL` |

**Return value**: `OCAS_OK` (0) on success.

**`OcasTensorContraction` fields**:
| Field | Type | Description |
|---|---|---|
| `kind` | `int` | `0` = product (free indices remain), `1` = scalar (fully contracted) |
| `tensors` | `OcasTensor**` | Valid when `kind == 0`; array of tensor handles |
| `n_tensors` | `size_t` | Length of the `tensors` array |
| `scalar_str` | `char*` | Valid when `kind == 1`; string of the scalar result |

**Memory**:
- `kind == 0`: free each `tensors[i]` with `ocas_tensor_free()`, and the array with `ocas_tensor_contraction_free()`.
- `kind == 1`: free `scalar_str` with `ocas_string_free()`, and the struct with `ocas_tensor_contraction_free()`.

**Example**:
```c
OcasTensor *A = ocas_tensor_create("A", "i,upper;j,lower", NULL, NULL);
OcasTensor *B = ocas_tensor_create("B", "j,upper;k,lower", NULL, NULL);
OcasTensorContraction result = {0};
if (ocas_tensor_contract(A, B, &result, NULL) == OCAS_OK) {
    if (result.kind == 0) {
        // free indices remain
        for (size_t i = 0; i < result.n_tensors; i++) {
            char *s = ocas_tensor_to_string(result.tensors[i], NULL);
            printf("result tensor: %s\n", s);
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

**Signature**:
```c
void ocas_tensor_contraction_free(OcasTensorContraction *c);
```

**Description**: Frees the `tensors` array and `scalar_str` of a contraction result. **Does not free** the individual tensor handles. `NULL`-safe.

---

### `ocas_tensor_canonicalize`

**Signature**:
```c
char *ocas_tensor_canonicalize(const char *expr_str, const char *specs_str,
                                const char *groups_str, int *err);
```

**Description**: Canonicalizes a tensor expression via graph isomorphism (0.22.0).

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `expr_str` | `const char*` | Tensor product string, e.g. `"T(i,j)*U(j,k)"` |
| `specs_str` | `const char*` | Comma-separated `name:sym` pairs, e.g. `"T:none,U:none"` |
| `groups_str` | `const char*` | May be `NULL`; or comma-separated `label:group` pairs |
| `err` | `int*` | May be `NULL` |

**Return value**: A heap-allocated canonical form string. Free with `ocas_string_free()`.

**Example**:
```c
char *canon = ocas_tensor_canonicalize(
    "T(i,j)*T(j,i)", "T:symmetric", NULL, NULL);
printf("canonical form: %s\n", canon);
ocas_string_free(canon);
```

---

### `ocas_young_project`

**Signature**:
```c
char *ocas_young_project(const char *expr_str, const char *tableau_str, int *err);
```

**Description**: Applies a Young projector to a tensor expression (0.22.0).

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `expr_str` | `const char*` | Tensor expression, e.g. `"f(a,b)"` |
| `tableau_str` | `const char*` | Comma-separated row lengths, e.g. `"2,1"` denotes □□/□ |
| `err` | `int*` | May be `NULL` |

**Return value**: A heap-allocated expanded expression string. Free with `ocas_string_free()`.

**Example**:
```c
// fully antisymmetric projection ((1,1,1) Young tableau, 3 indices)
char *proj = ocas_young_project("f(a,b,c)", "1,1,1", NULL);
printf("projection result: %s\n", proj);
ocas_string_free(proj);
```

---

### `ocas_tensor_refresh_dummies`

**Signature**:
```c
char *ocas_tensor_refresh_dummies(const char *expr_str, const char *specs_str, int *err);
```

**Description**: Renames the dummy indices (labels occurring exactly twice) of a tensor expression to `d0`, `d1`, ... to avoid conflicts (0.22.0).

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `expr_str` | `const char*` | Tensor expression |
| `specs_str` | `const char*` | `name:sym` pairs |
| `err` | `int*` | May be `NULL` |

**Return value**: A heap-allocated string. Free with `ocas_string_free()`.

---

### C++ wrapper: `ocas::tensor::*`

```cpp
std::string canon = ocas::tensor::canonicalize("T(i,j)*T(j,i)", "T:symmetric");
std::string proj = ocas::tensor::young_project("f(a,b,c)", "2,1");
std::string fresh = ocas::tensor::refresh_dummies("T(i,j)*U(j,i)", "T:none,U:none");
```

All wrappers throw `ocas::Error` on error.

---

## 11. Dual Numbers (14 functions)

Hyper-dual numbers (forward automatic differentiation), supporting only `Rational` coefficients. Polynomial/rational arithmetic.

### Coefficient string format

Rationals are passed as `"num"` (denominator 1) or `"num/den"`. Returned strings use the same format.

### `ocas_dual_shape_new`

**Signature**:
```c
OcasDualShape *ocas_dual_shape_new(size_t n_vars, int *err);
```

**Description**: Creates a first-order shape tracking derivatives with respect to `n_vars` variables.

**Return value**: `OcasDualShape*` handle or `NULL`.

**Example**:
```c
OcasDualShape *shape = ocas_dual_shape_new(2, NULL);  // 2 variables
printf("n vars: %zu, n components: %zu\n",
       ocas_dual_shape_n_vars(shape),
       ocas_dual_shape_n_components(shape));
// Output: n vars: 2, n components: 3
ocas_dual_shape_free(shape);
```

---

### `ocas_dual_shape_free`

**Signature**:
```c
void ocas_dual_shape_free(OcasDualShape *s);
```

**Description**: Frees a shape handle. `NULL`-safe.

---

### `ocas_dual_shape_n_vars`

**Signature**:
```c
size_t ocas_dual_shape_n_vars(const OcasDualShape *s);
```

**Description**: Returns the number of differentiation variables. `NULL` returns `0`.

---

### `ocas_dual_shape_n_components`

**Signature**:
```c
size_t ocas_dual_shape_n_components(const OcasDualShape *s);
```

**Description**: Returns the total number of components (value + each partial derivative). `NULL` returns `0`.

---

### `ocas_dual_variable`

**Signature**:
```c
OcasHyperDual *ocas_dual_variable(const OcasDualShape *shape, size_t i,
                                  const char *coeff, int *err);
```

**Description**: Creates an independent variable $x_i = \text{coeff}$ (derivative with respect to variable $i$ is 1).

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `shape` | `const OcasDualShape*` | Non-NULL shape handle |
| `i` | `size_t` | Variable index ($0 \le i < \text{n\_vars}$) |
| `coeff` | `const char*` | Coefficient string `"num"` or `"num/den"` |
| `err` | `int*` | May be `NULL` |

**Return value**: `OcasHyperDual*` handle or `NULL`.

---

### `ocas_dual_constant`

**Signature**:
```c
OcasHyperDual *ocas_dual_constant(const OcasDualShape *shape,
                                  const char *coeff, int *err);
```

**Description**: Creates a constant dual number (all derivatives zero).

---

### `ocas_hyperdual_free`

**Signature**:
```c
void ocas_hyperdual_free(OcasHyperDual *d);
```

**Description**: Frees a hyper-dual number handle. `NULL`-safe.

---

### `ocas_dual_value`

**Signature**:
```c
char *ocas_dual_value(const OcasHyperDual *d, int *err);
```

**Description**: Returns a heap-allocated string of the scalar value component. Free with `ocas_string_free()`.

---

### `ocas_dual_deriv`

**Signature**:
```c
char *ocas_dual_deriv(const OcasHyperDual *d, size_t i, int *err);
```

**Description**: Returns a heap-allocated string of the partial derivative with respect to variable $i$. Returns `NULL` if the shape has no first-order component for $i$.

**Example**:
```c
OcasDualShape *shape = ocas_dual_shape_new(2, NULL);
OcasHyperDual *x = ocas_dual_variable(shape, 0, "3", NULL);   // x₀ = 3
OcasHyperDual *y = ocas_dual_variable(shape, 1, "5", NULL);   // x₁ = 5
OcasHyperDual *prod = ocas_dual_mul(x, y, NULL);               // f = x₀·x₁

char *val = ocas_dual_value(prod, NULL);
char *df_dx = ocas_dual_deriv(prod, 0, NULL);  // ∂f/∂x₀ = x₁ = 5
char *df_dy = ocas_dual_deriv(prod, 1, NULL);  // ∂f/∂x₁ = x₀ = 3

printf("f = %s, ∂f/∂x₀ = %s, ∂f/∂x₁ = %s\n", val, df_dx, df_dy);
// Output: f = 15, ∂f/∂x₀ = 5, ∂f/∂x₁ = 3

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

**Signature**:
```c
OcasHyperDual *ocas_dual_add(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**Description**: Computes $a + b$. Both operands must share the same shape.

**Return value**: A new handle or `NULL`.

---

### `ocas_dual_sub`

**Signature**:
```c
OcasHyperDual *ocas_dual_sub(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**Description**: Computes $a - b$.

---

### `ocas_dual_mul`

**Signature**:
```c
OcasHyperDual *ocas_dual_mul(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**Description**: Computes $a \times b$. Derivatives are propagated automatically via the product rule.

---

### `ocas_dual_div`

**Signature**:
```c
OcasHyperDual *ocas_dual_div(const OcasHyperDual *a, const OcasHyperDual *b, int *err);
```

**Description**: Computes $a / b$. If the value component of $b$ is zero, sets `OCAS_ERROR_DIVISION_BY_ZERO` and returns `NULL`.

---

### `ocas_dual_neg`

**Signature**:
```c
OcasHyperDual *ocas_dual_neg(const OcasHyperDual *a, int *err);
```

**Description**: Computes $-a$.

---

## 12. Gröbner (8 functions)

Gröbner basis computation and polynomial system solving.

> **⚠️ Non-standard conventions**
>
> - These functions are **not included in the `ocas.h` header**; declare the prototypes yourself.
> - Error codes use `-1` instead of the `OCAS_*` constants.
> - The `OcasGroebnerBasis` and `OcasSystemSolution` handles are **not marked** `#[repr(C)]`.

### Polynomial data format

Polynomials are passed as arrays of coefficients; each polynomial is specified by:
- `n_vars` — number of variables
- `n_terms` — number of terms
- `exponents` — flattened exponent matrix (`n_terms × n_vars`, row-major)
- `coeff_nums` / `coeff_dens` — numerator/denominator arrays of the rational coefficients

### `ocas_groebner_basis`

**Signature**:
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

**Description**: Computes the Gröbner basis from an array of polynomial data.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `n_polys` | `size_t` | Number of polynomials |
| `n_vars_array` | `const size_t*` | Number of variables for each polynomial |
| `n_terms_array` | `const size_t*` | Number of terms for each polynomial |
| `exponents` | `const size_t*` | Flattened exponent matrix |
| `coeff_nums` | `const int64_t*` | Numerators of the rational coefficients |
| `coeff_dens` | `const int64_t*` | Denominators of the rational coefficients |
| `algorithm` | `int32_t` | Algorithm selection: `0` = Auto, `1` = F4, `2` = F5, `3` = Buchberger |
| `err` | `int32_t*` | Error code (`-1` means failure). **Note**: unlike the rest of the C API, this pointer must NOT be `NULL` — the function returns `NULL` without computing when it is |

**Return value**: `OcasGroebnerBasis*` handle or `NULL`.

**Memory**: The caller frees with `ocas_groebner_basis_free()`.

**Example**:
```c
// Compute the Gröbner basis of {x² + y - 1, x + y² - 1}
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
    printf("the basis has %zu elements\n", ocas_groebner_basis_len(gb));
    printf("zero-dimensional? %s\n", ocas_is_zero_dimensional(gb) ? "yes" : "no");
    ocas_groebner_basis_free(gb);
}
```

---

### `ocas_groebner_basis_free`

**Signature**:
```c
void ocas_groebner_basis_free(OcasGroebnerBasis *gb);
```

**Description**: Frees a Gröbner basis handle. `NULL`-safe.

---

### `ocas_groebner_basis_len`

**Signature**:
```c
size_t ocas_groebner_basis_len(const OcasGroebnerBasis *gb);
```

**Description**: Returns the number of elements in the Gröbner basis. `NULL` returns `0`.

---

### `ocas_is_zero_dimensional`

**Signature**:
```c
bool ocas_is_zero_dimensional(const OcasGroebnerBasis *gb);
```

**Description**: Checks whether the ideal is zero-dimensional. `NULL` returns `false`.

---

### `ocas_solve_polynomial_system`

**Signature**:
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

**Description**: Solves a system of polynomial equations. Returns a handle or `NULL`.

**Parameters**: Same as `ocas_groebner_basis`.

**Memory**: The caller frees with `ocas_system_solution_free()`.

---

### `ocas_system_solution_count`

**Signature**:
```c
size_t ocas_system_solution_count(const OcasSystemSolution *sol);
```

**Description**: Returns the number of solutions. Returns `0` for positive-dimensional or empty solution sets. `NULL` returns `0`.

---

### `ocas_system_solution_value`

**Signature**:
```c
double ocas_system_solution_value(const OcasSystemSolution *sol,
                                  size_t sol_idx, size_t var_idx);
```

**Description**: Gets the value of a specific variable of a specific solution.

**Parameters**:
| Parameter | Type | Description |
|---|---|---|
| `sol` | `const OcasSystemSolution*` | Non-NULL solution handle |
| `sol_idx` | `size_t` | Index of the solution |
| `var_idx` | `size_t` | Index of the variable |

**Return value**: An `f64` value. Returns `0.0` on out-of-range or error.

**Example**:
```c
OcasSystemSolution *sol = ocas_solve_polynomial_system(
    2, n_vars, n_terms, exps, nums, dens, 0, &err);
if (sol) {
    size_t count = ocas_system_solution_count(sol);
    printf("found %zu solutions\n", count);
    for (size_t i = 0; i < count; i++) {
        printf("solution %zu: x = %g, y = %g\n", i,
               ocas_system_solution_value(sol, i, 0),
               ocas_system_solution_value(sol, i, 1));
    }
    ocas_system_solution_free(sol);
}
```

---

### `ocas_system_solution_free`

**Signature**:
```c
void ocas_system_solution_free(OcasSystemSolution *sol);
```

**Description**: Frees a solution handle. `NULL`-safe.

---

## Complete Examples

### C: Symbolic differentiation and evaluation

```c
#include <ocas.h>
#include <stdio.h>
#include <assert.h>

int main(void) {
    int err = 0;

    // Parse the expression
    OcasExpr *f = ocas_expr_parse("x^3 - 6*x^2 + 11*x - 6", &err);
    assert(f != NULL);

    // Differentiate
    OcasExpr *df = ocas_expr_diff(f, "x", &err);
    assert(df != NULL);

    // Simplify
    OcasExpr *simplified = ocas_expr_simplify(df, &err);
    assert(simplified != NULL);

    char *s = ocas_expr_to_string(simplified, &err);
    printf("f'(x) = %s\n", s);  // Output: f'(x) = 11 + (-12*x) + (3*(x^2))

    // Cleanup
    ocas_string_free(s);
    ocas_expr_free(simplified);
    ocas_expr_free(df);
    ocas_expr_free(f);

    return 0;
}
```

Compile:
```bash
gcc example.c -locas -o example
```

### C: Automatic differentiation (dual numbers)

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // f(x,y) = x² + 2xy + y², evaluated at (3, 5) with its partial derivatives
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

> **Note**: the dual-number domain has **no** C++ RAII wrappers in `ocas.hpp`; use the C API above directly.

### C: Number theory

```c
#include <ocas.h>
#include <stdio.h>

int main(void) {
    int err = 0;

    // Prime factorization
    char *f = ocas_ntheory_factorint("360", &err);
    printf("360 = %s\n", f);  // 2:3,3:2,5:1
    ocas_string_free(f);

    // Primality test
    printf("is 997 prime? %s\n", ocas_ntheory_isprime("997", &err) ? "yes" : "no");

    // Chinese remainder theorem
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

## See Also

- [Rust API reference](./rust.md) — overview of the Rust public API
- [Python API reference](./python.md) — Python binding documentation
- [Architecture](../architecture.md) — overview of the project architecture
