# Getting Started

This chapter gets you up and running with oCAS in about 5 minutes. oCAS 0.24 supports three language interfaces: **Rust**, **Python**, and **C/C++**.

> **System requirements**
>
> | Language | Minimum version | Notes |
> |---|---|---|
> | Rust | 1.97+ (edition 2024) | Install via [rustup](https://rustup.rs/) recommended |
> | Python | 3.10+ | Use `pip install` or build from source |
> | C/C++ | C11 / C++17 | Building the C library requires `cargo` |

---

## Rust

### Installation

Add oCAS to your `Cargo.toml`:

```toml
[dependencies]
ocas = "0.24"
```

To enable optional backends, turn on features:

```toml
[dependencies]
ocas = { version = "0.24", features = ["gmp", "jit"] }
```

Commonly used features:

| Feature | Description |
|---|---|
| `gmp` | Use GMP arbitrary-precision arithmetic (instead of pure-Rust `num-bigint`) |
| `mpfr` | Use MPFR ball arithmetic (real interval arithmetic) |
| `jit` | Cranelift JIT-compiled evaluator |
| `simd` | SIMD-vectorized batch evaluation |
| `egg` | e-graph equality saturation engine |
| `ntt` | NTT fast polynomial multiplication (𝔽_p) |
| `mimalloc` | Use the mimalloc allocator |
| `system-libs` | Link against system-installed GMP/MPFR (required on Windows/MinGW) |

See the [Rust API reference](./api/rust.md) for the complete feature list.

### Verify the Installation

Create a minimal project to confirm everything works:

```bash
cargo new ocas-test && cd ocas-test
```

Add the dependency to `Cargo.toml`:

```toml
[dependencies]
ocas = "0.24"
```

Write the following to `src/main.rs`:

```rust
use ocas::prelude::*;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);
    let expr = parse(&ctx, "x^2 + 2*x + 1").unwrap();
    let d = diff(&ctx, expr, Symbol::new("x"));
    println!("d/dx(x² + 2x + 1) = {}", d);
    // Output: d/dx(x² + 2x + 1) = 2 + 2*x
}
```

Run it:

```bash
cargo run
```

If it prints `2 + 2*x`, the installation succeeded.

### Basic Usage

#### Symbolic Differentiation and Simplification

```rust
use ocas::prelude::*;
use ocas::ocas_rewrite::rules::default_rules;

fn main() {
    let arena = Arena::new();
    let ctx = AtomArena::new(&arena);

    // Parse the expression
    let expr = parse(&ctx, "x^2").unwrap();

    // Differentiate
    let d = diff(&ctx, expr, Symbol::new("x"));
    println!("derivative: {}", d);
    // Output: 2*x

    // Simplify (default rule set: identity removal, constant folding)
    let e2 = parse(&ctx, "x + 0 + y*0 + z*1").unwrap();
    let rules = default_rules(&ctx, &());
    let s = simplify(&ctx, e2, &rules, 100);
    println!("simplified: {}", s);
    // Output: x + z
}
```

#### Polynomial Arithmetic

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
    // DenseUnivariatePolynomial does not implement Display; assert coefficients
    assert_eq!(quot.coeffs(), &[Rational::new(1, 1), Rational::new(1, 1)]); // 1 + x
    assert_eq!(rem.coeffs(), &[]); // remainder 0
}
```

#### Solving Linear Systems

```rust
use ocas_calc::solve::solve_linear_rational;

fn main() {
    // Solve 2x + y = 5, x - y = 1
    let a = vec![vec![2, 1], vec![1, -1]];
    let b = vec![5, 1];

    let solution = solve_linear_rational(&a, &b).unwrap();
    println!("x = {}/{}, y = {}/{}",
        solution[0].0, solution[0].1,
        solution[1].0, solution[1].1);
    // Output: x = 2/1, y = 1/1
}
```

#### Number Theory

```rust
use ocas_domain::number_theory::{
    is_prime, factor::factor_integer, functions::euler_phi,
};

fn main() {
    println!("is_prime(97) = {}", is_prime(&97.into()));
    // Output: is_prime(97) = true

    let factors = factor_integer(&360.into());
    println!("360 = {}", factors.iter()
        .map(|(p, e)| format!("{}^{}", p, e))
        .collect::<Vec<_>>()
        .join(" × "));
    // Output: 360 = 2^3 × 3^2 × 5^1

    println!("φ(360) = {}", euler_phi(&360.into()));
    // Output: φ(360) = 96
}
```

For more Rust APIs, see:
- [Expression system](./api/rust-expressions.md) — Arena, Symbol, Atom, pattern matching
- [Coefficient domains](./api/rust-domains.md) — Integer, Rational, FiniteField, etc.
- [Polynomials](./api/rust-polynomials.md) — univariate/multivariate polynomials, monomial orders
- [Matrices](./api/rust-matrix.md) — determinant, solving, multiplication
- [Calculus](./api/rust-calculus.md) — diff, integrate, taylor
- [Solvers](./api/rust-solvers.md) — linear equations, Diophantine, ODEs
- [Number theory](./api/rust-ntheory.md) — primality, factorization, discrete logarithms

---

## Python

### Installation

```bash
pip install ocas
```

> Requires Python ≥ 3.10 (uses the PyO3 abi3 stable ABI).

### Verify the Installation

```bash
python -c "import ocas; e = ocas.Expression('x^2'); print(e.diff('x'))"
```

If it prints `2*x`, the installation succeeded. To verify in an interactive Python session:

```python
>>> import ocas
>>> e = ocas.Expression("x^2 + 2*x + 1")
>>> e.diff("x")
2 + 2*x
>>> e.simplify()
1 + 2*x + x^2
```

### Basic Usage

#### Symbolic Differentiation and Simplification

```python
import ocas

e = ocas.Expression("x^2")

# Differentiate
d = e.diff("x")
print(f"derivative: {d}")
# Output: derivative: 2*x

# Simplify (default rule set: identity removal, constant folding)
e2 = ocas.Expression("x + 0 + y*0 + z*1")
s = e2.simplify()
print(f"simplified: {s}")
# Output: simplified: x + z
```

#### Polynomial Arithmetic

```python
import ocas

# 1 + 2x + x²
p = ocas.Polynomial([1, 2, 1])
print(f"degree: {p.degree()}")         # degree: 2
print(f"p(2) = {p.eval(2)}")         # p(2) = 9

# Polynomial GCD
q = ocas.Polynomial([1, 1])          # 1 + x
g = p.gcd(q)
print(f"gcd = {g}")                   # gcd = 1 + x

# Factorization
factors = p.factor()
print(f"factors: {factors}")         # factors: [(1 + x, 2)]
```

#### Finite Field Arithmetic

```python
import ocas

gf5 = ocas.FiniteField(5)
q = ocas.Polynomial([1, 2, 1], domain=gf5)
print(q.eval(3))  # (1 + 6 + 9) mod 5 = 16 mod 5 = 1
```

#### Matrix Operations

```python
import ocas

m = ocas.Matrix([[1, 2], [3, 4]])
print(f"determinant: {m.determinant()}")       # determinant: -2
print(f"transpose: {(m.transpose()).rows()}")   # transpose: [[1, 3], [2, 4]]
print(f"m × m: {(m @ m).rows()}")          # m × m: [[7, 10], [15, 22]]
```

#### Number Theory

```python
import ocas

print(ocas.isprime(97))               # True
print(ocas.factorint(360))             # [("2", 3), ("3", 2), ("5", 1)]
print(ocas.totient(360))               # 96
print(ocas.nextprime(100))             # 101
print(ocas.mobius(30))                 # -1
print(ocas.divisor_count(360))         # 24
```

#### ODE Solving

```python
import ocas

e = ocas.Expression("Derivative(y(x), x) + y(x)")  # y' + y = 0

# Classify the available solving methods
methods = ocas.classify_ode(e, "y", "x")
print(methods)

# Solve y' + y = 0
sol = ocas.dsolve(e, "y", "x")
print(sol)  # y = C1*exp(-x)
```

For more Python APIs, see the [Python API reference](./api/python.md).

---

## C/C++

### Building the C Library

oCAS's C library is built with `cargo`. Clone the source repository and compile:

```bash
git clone https://github.com/charleshzh/ocas.git
cd ocas
cargo build -p ocas-c --release
```

The build artifacts are placed in `target/release/`:

| File | Description |
|---|---|
| `libocas_c.so` (Linux) / `ocas_c.dll` (Windows) / `libocas_c.dylib` (macOS) | Dynamic library |
| `libocas_c.a` (Linux/macOS) / `ocas_c.lib` (Windows) | Static library |

The header file `ocas.h` lives in the `ocas-c/include/` directory. C++ users can also use `ocas.hpp` in the same directory (RAII wrappers covering the expression, number theory, and tensor domains).

> **Windows users**: the MSYS2 MINGW64 environment is recommended. See [Building on Windows](./build-windows.md).

### Verify the Installation

Create a minimal test file `test_ocas.c`:

```c
#include <stdio.h>
#include "ocas.h"

int main(void) {
    struct ocas_OcasExpr *e = ocas_expr_parse("x^2 + 2*x + 1", NULL);
    if (!e) {
        fprintf(stderr, "parse failed\n");
        return 1;
    }

    struct ocas_OcasExpr *d = ocas_expr_diff(e, "x", NULL);
    char *s = ocas_expr_to_string(d, NULL);
    printf("d/dx(x² + 2x + 1) = %s\n", s);
    /* Output: d/dx(x² + 2x + 1) = 2 + 2*x */

    ocas_string_free(s);
    ocas_expr_free(d);
    ocas_expr_free(e);
    return 0;
}
```

Compile and run (Linux example; adjust the paths to match your environment):

```bash
# assuming oCAS is built at /path/to/ocas
gcc test_ocas.c -I/path/to/ocas/ocas-c/include \
    -L/path/to/ocas/target/release -locas_c -lm -ldl -lpthread -o test_ocas
LD_LIBRARY_PATH=/path/to/ocas/target/release ./test_ocas
```

If it prints `2 + 2*x`, the installation succeeded.

### Basic Usage

#### Symbolic Differentiation

```c
#include <stdio.h>
#include "ocas.h"

int main(void) {
    struct ocas_OcasExpr *e = ocas_expr_parse("x^2", NULL);
    struct ocas_OcasExpr *d = ocas_expr_diff(e, "x", NULL);

    char *s = ocas_expr_to_string(d, NULL);
    printf("derivative: %s\n", s);
    /* Output: derivative: 2*x */

    ocas_string_free(s);
    ocas_expr_free(d);
    ocas_expr_free(e);
    return 0;
}
```

#### Number Theory

```c
#include <stdio.h>
#include "ocas.h"

int main(void) {
    char *result;

    /* Primality test */
    int prime = ocas_ntheory_isprime("97", NULL);
    printf("is_prime(97) = %s\n", prime ? "true" : "false");
    /* Output: is_prime(97) = true */

    /* Integer factorization */
    result = ocas_ntheory_factorint("360", NULL);
    printf("360 = %s\n", result);
    /* Output: 360 = 2:3,3:2,5:1 */
    ocas_string_free(result);

    /* Euler φ function */
    result = ocas_ntheory_totient("360", NULL);
    printf("φ(360) = %s\n", result);
    /* Output: φ(360) = 96 */
    ocas_string_free(result);

    return 0;
}
```

#### C++ RAII Wrappers

If you use C++, the RAII wrappers in `ocas.hpp` save you from manual memory management:

```cpp
#include <iostream>
#include "ocas.h"
#include "ocas.hpp"

int main() {
    ocas::Expression e("x^2 + 2*x + 1");
    auto d = e.diff("x");
    std::cout << "d/dx = " << d.to_string() << std::endl;
    // Output: d/dx = 2 + 2*x

    // Number theory (no manual free)
    std::cout << "is_prime(97) = " << ocas::ntheory::isprime("97") << std::endl;
    std::cout << "360 = " << ocas::ntheory::factorint("360") << std::endl;
    std::cout << "φ(360) = " << ocas::ntheory::totient("360") << std::endl;

    return 0;
}
```

For more C/C++ APIs, see the [C/C++ API reference](./api/c.md).

---

## FAQ

### Rust

**Q: Compilation is slow, what can I do?**

oCAS contains a lot of generic code, so the first build can be slow. Suggestions:
- During development use `cargo build` (dev profile, dependencies already use `opt-level = 3`)
- For releases use `cargo build --release` (enables LTO)
- Use `sccache` or `cargo-nextest` to speed up incremental builds

**Q: The `gmp` feature fails because the GMP library is not found?**

On Linux/macOS, make sure the GMP development package is installed:
```bash
# Ubuntu/Debian
sudo apt install libgmp-dev

# macOS (Homebrew)
brew install gmp
```

On Windows (MSYS2 MINGW64):
```bash
pacman -S mingw-w64-x86_64-gmp
```

Or use the `system-libs` feature to link against system-installed libraries. See [Building on Windows](./build-windows.md).

**Q: The `flint` feature is unavailable on Windows?**

The upstream `flint3-sys` crate only supports POSIX platforms, so the `flint` feature is not supported on Windows for now.

### Python

**Q: `pip install ocas` fails?**

Check that Python ≥ 3.10 is used:
```bash
python --version
```

To build from source:
```bash
pip install maturin
maturin develop --release -m ocas-py/Cargo.toml
```

**Q: Why do some functions in the Python bindings have a `py_` prefix?**

Some Gröbner functions (e.g. `ocas.py_groebner_basis`) keep the `py_` prefix because the PyO3 bindings have not yet added `#[pyo3(name=...)]` annotations. Just call the functions with the prefixed names.

### C/C++

**Q: Linking fails with undefined reference?**

Make sure the link order is correct — `libocas_c` goes after your source file:
```bash
gcc test.c -locas_c -lm -ldl -lpthread
```

And make sure `LD_LIBRARY_PATH` (Linux) or `PATH` (Windows) contains the library directory.

**Q: Some functions are missing from `ocas.h`?**

The `ocas.h` shipped in the repository may not be up to date. cbindgen generates a fresh header under `target/release/build/ocas-c-*/out/` at build time. Use that file instead.

**Q: What functionality do the C++ wrappers in `ocas.hpp` cover?**

`ocas.hpp` currently covers three domains: **expressions** (`ocas::Expression`), **number theory** (`ocas::ntheory::*`), and **tensors** (`ocas::tensor::*`). For the remaining domains (polynomials, Gröbner, ODEs, etc.) use the C API directly.

### General

**Q: Numeric results are imprecise?**

By default oCAS uses pure-Rust `num-bigint` rational arithmetic, which is exact. If you need high-precision floats (interval arithmetic), enable the `mpfr` feature:
```toml
ocas = { version = "0.24", features = ["mpfr"] }
```

**Q: How do I check the oCAS version?**

```bash
cargo metadata --format-version=1 | jq '.packages[] | select(.name=="ocas") | .version'
```

---

## Next Steps

- [Rust API reference](./api/rust.md) — complete documentation of Rust modules, types, and functions
- [Python API reference](./api/python.md) — detailed documentation of all Python classes and functions
- [C/C++ API reference](./api/c.md) — C exported functions, C++ RAII wrappers, memory management
- [Architecture](./architecture.md) — overall design and module relationships of oCAS
- [Gröbner basis implementation](./algorithms/groebner.md) — Buchberger / F4 / F5 algorithms in detail
- [ODE solving](./algorithms/ode-solving.md) — classification and solution methods for ordinary differential equations
