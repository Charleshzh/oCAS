# ODE Solving Algorithms

oCAS's ODE module (`ocas-calc/src/ode/`) implements a symbolic solver for ordinary differential equations. This chapter details the classification criteria for each ODE type, the implementation details of the solving algorithms, and the design choices. For the mathematical derivations, see [ODE Solving Theory](../math/ode-theory.md).

---

## Overall Architecture

The core entry point of ODE solving is `dsolve`, which performs three steps:

1. **Normalization** (`normalize_ode`): algebraically simplifies and standardizes the input equation
2. **Classification** (`classify_ode`): returns all applicable solving methods in priority order
3. **Dispatch**: tries the methods in turn and returns the first successful result

```
input equation ─→ normalize ─→ classify ─→ dispatch ─→ ODESolution
                              │
          ┌───────────────────┼───────────────────┐
          │                   │                   │
     first-order ODE    second-order linear ODE   series solutions (fallback)
    ┌────┴────┐          ┌────┴────┐
 separable  linear   const-coeff  Cauchy-Euler
 Bernoulli  exact    reduction    variation of
 homogeneous          of order    parameters
```

### Result Types

The `ODESolution` enum represents the outcome:

| Variant | Meaning | Typical source |
|---|---|---|
| `Explicit(expr)` | explicit solution $y = f(x)$ | linear, constant coefficient, Laplace |
| `Implicit(expr)` | implicit solution $F(x, y) = C$ | separable, exact, homogeneous |
| `Parametric(x(t), y(t))` | parametric solution $(x(t), y(t))$ | reserved in the enum; no solver produces it yet |
| `Series(expr, n)` | truncated series solution | power series, Frobenius |
| `System(&[Atom])` | system solution components $(y_1, y_2)$ | $2\times2$ systems |
| `Unsolved(ode)` | could not be solved symbolically | when all methods fail |

---

## The Classification Engine

The classifier (`classify.rs`) performs a structural analysis of the ODE and returns the list of candidate methods in priority order. It does not solve — only pattern recognition.

### Detection Pipeline

**First-order ODEs** (`order == 1`):

1. Check `is_first_order_linear`: $y$ and $y'$ appear only to the first power → `LinearFirst`
2. Check `is_bernoulli`: some $y^n$ ($n \geq 2$) plus a linear $y$ term → `Bernoulli`
3. Check `is_separable`: additive terms without $y$ plus additive terms with $y$ → `Separable`
4. Check `is_exact`: in $M + Ny' = 0$, $\partial M/\partial y = \partial N/\partial x$ → `Exact`
5. Check `is_homogeneous`: all terms have the same total $(x, y)$ degree → `Homogeneous`

**Second-order and higher linear ODEs**:

6. `is_constant_coeff_linear`: the coefficients of $y, y', y''$ do not contain $x$ → `LinearConstantCoeff`
7. `is_cauchy_euler`: the terms are of the form $c_k x^k y^{(k)}$ → `CauchyEuler`
8. Second-order linear → `ReductionOfOrder` (tries simple candidate solutions)

**Global fallback**:

9. Linear ODEs → `PowerSeries` (power series or Frobenius)

### Exactness Detection

The exact-equation test `is_exact` works as follows:

1. Split the equation into $M(x, y) + N(x, y)\,y' = 0$ with `split_mn`
2. Replace $y(x)$ by the bare symbol $y$ so that $M, N$ become bivariate expressions
3. Compute $\partial M/\partial y$ and $\partial N/\partial x$
4. First compare the normalized strings (`to_string()` after `normalize`)
5. If they differ, compute the difference and check with the rule-based simplifier whether it is zero

```rust
let dm_norm = normalize(ctx, dm_dy);
let dn_norm = normalize(ctx, dn_dx);
if dm_norm.to_string() == dn_norm.to_string() {
    return true;
}
// fallback: simplify the difference
let difference = simplify(ctx, ctx.add(&[dm_dy, ctx.mul(&[ctx.num(-1), dn_dx])]), &rules, 20);
matches!(normalize(ctx, difference).node(), AtomNode::Num(0))
```

This double-check strategy compensates for the simplifier's inability to collect like terms.

### Homogeneity Detection

`is_homogeneous` checks whether the equation is homogeneous in $(x, y)$ (zero-th degree homogeneity). In the implementation, `all_terms_homogeneous_degree` walks all additive terms, computes the total degree of each term in $x$ and `func` (`term_degree`), and asserts that all terms have the same, positive degree.

For $\text{Derivative}(y(x), x)$, the degree is counted as 1 in $y$ (without additional powers of $x$).

---

## First-Order ODE Solvers

The first-order solvers are implemented in `first_order.rs`; there are five methods.

### Separable Equations

**Criterion**: the equation can be written as $g(y)\,y' = f(x)$.

**Algorithm**:

1. Split the equation into $y$-terms and non-$y$ ($x$) terms with `separate_by_func`
2. Replace $x$ by $y$ in the $y$-terms (`substitute_var`)
3. Integrate both sides: $\int g(y)\,dy$ and $\int f(x)\,dx$
4. Return the implicit solution $\int g(y)\,dy - \int f(x)\,dx = C$

**Design choice**: `Implicit` is returned rather than `Explicit` because in general solving for $y$ is required, which is undecidable in symbolic computation. The integration calls `crate::integral::integrate`, which launches the full integration pipeline (fast table → rational function → Risch → trigonometric rewrite → special functions → heuristic techniques → unevaluated form).

### First-Order Linear Equations

**Criterion**: $y' + p(x)\,y = q(x)$, where $y$ and $y'$ appear only to the first power.

**Algorithm** (integrating factor method):

1. Extract $p(x)$ and $q(x)$ with `extract_linear_coeffs`
2. Compute the integrating factor $\mu(x) = e^{\int p\,dx}$
3. General solution:
$$y = e^{-\int p\,dx}\left(\int e^{\int p\,dx}\,q(x)\,dx + C_1\right)$$

```rust
let p_integral = integrate(ctx, p, var);
let mu = ctx.fun("exp", &[p_integral]);
let mu_q = ctx.mul(&[mu, q]);
let int_mu_q = integrate(ctx, mu_q, var);
let neg_p_integral = ctx.mul(&[ctx.num(-1), p_integral]);
let exp_neg_p = ctx.fun("exp", &[neg_p_integral]);
let particular = ctx.mul(&[exp_neg_p, int_mu_q]);
let c1 = ctx.var("C1");
let homogeneous = ctx.mul(&[c1, exp_neg_p]);
```

`extract_linear_coeffs` decomposes the equation into additive terms and classifies them as $y'$-terms, $y$-terms, and free terms. For the $y$-terms it extracts the coefficient $p(x)$ by filtering out the `func` factor.

### Bernoulli Equations

**Criterion**: $y' + p(x)\,y = q(x)\,y^n$ ($n \neq 0, 1$).

**Algorithm**:

1. Find $n$ in $y^n$ recursively with `find_bernoulli_power`
2. Substitute $v = y^{1-n}$; the equation becomes $v' + (1-n)\,p(x)\,v = (1-n)\,q(x)$
3. Solve for $v$ with the linear-equation method
4. Substitute back $y = v^{1/(1-n)}$

`find_power_inner` walks the AST recursively and returns $n$ when it encounters `Pow(base, Num(n))` whose `base` contains `func`.

**Edge cases**: $n = 0$ or $n = 1$ degenerates to a linear equation; `None` is returned so that `LinearFirst` handles it.

### Exact Equations

**Criterion**: $M(x,y) + N(x,y)\,y' = 0$ with $\partial M/\partial y = \partial N/\partial x$.

**Algorithm**:

1. Split into $M$ and $N$ with `split_mn` ($y'$ must appear linearly, otherwise `None` is returned)
2. Check exactness (`partials_equal`)
3. If not exact, try an integrating factor (`find_integrating_factor`)
4. Integrate $M$ with respect to $x$ to obtain part of the potential $F$
5. Compute the correction term $g(y) = N - \partial F/\partial y$
6. If $g(y)$ does not contain $x$, integrate with respect to $y$ to obtain the full solution $F + g = C$

**Integrating factor strategies**:

- **Candidate 1**: $(M_y - N_x)/N$ contains only $x$ → $\mu(x) = \exp\!\int\!(M_y - N_x)/N\,dx$
- **Candidate 2**: $(N_x - M_y)/M$ contains only $y$ → $\mu(y) = \exp\!\int\!(N_x - M_y)/M\,dy$

The helper `exp_simplify` simplifies $\exp(k \cdot \log u) = u^k$, so that the integrating factor does not remain in unsimplified exponential form.

### Homogeneous Equations

**Criterion**: $y' = f(y/x)$ (all additive terms have the same total $(x,y)$ degree).

**Algorithm**:

1. Substitute $v = y/x$ (i.e. $y = vx$), obtaining $y' = v + xv'$
2. The equation becomes $v + xv' = f(v)$, i.e. $xv' = f(v) - v$ (separable in $v$ and $x$)
3. Split into $v$-terms and $x$-terms with `separate_by_var`
4. Integrate both sides, then substitute back $v = y/x$

Returns an `Implicit` solution.

---

## Second-Order Linear ODE Solvers

The second-order solvers are implemented in `second_order.rs`.

### Constant-Coefficient Equations

**Form**: $a\,y'' + b\,y' + c\,y = f(x)$, where $a, b, c$ are constants.

**Algorithm**:

1. Extract $a, b, c$ and $f(x)$ with `extract_second_order_coeffs`
2. Verify that $a, b, c$ do not contain $x$ (`contains_x`)
3. Compute the discriminant $\Delta = b^2 - 4ac$ of the characteristic equation $ar^2 + br + c = 0$
4. Construct the homogeneous solution $y_c$
5. Construct a particular solution $y_p$ (undetermined coefficients or variation of parameters)
6. Return $y = y_c + y_p$

#### Three Cases for the Homogeneous Solution

`constant_coeff_basis` constructs the fundamental solution set from the discriminant:

| Discriminant | Fundamental solution set |
|---|---|
| $\Delta > 0$ (two real roots $r_1, r_2$) | $e^{r_1 x},\; e^{r_2 x}$ |
| $\Delta = 0$ (repeated root $r$) | $e^{rx},\; x\,e^{rx}$ |
| $\Delta < 0$ (complex roots $\alpha \pm \beta i$) | $e^{\alpha x}\cos\beta x,\; e^{\alpha x}\sin\beta x$ |

In the implementation, the discriminant branches directly once it has been collapsed (`collect_terms`) to `Num(d)`. For complex roots, $\beta = \sqrt{-\Delta}/(2a)$: when $-\Delta$ is a perfect square, $\beta$ is computed exactly as a rational number; otherwise the symbolic radical is kept.

```rust
let sn = isqrt(num);
if sn * sn == num {
    // beta = sn/den — exact rational
    ctx.mul(&[ctx.num(sn), ctx.pow(ctx.num(den), ctx.num(-1))])
} else {
    // beta = sqrt(num)/den — symbolic form
    ctx.mul![
        ctx.pow(ctx.num(num), ctx.pow(ctx.num(2), ctx.num(-1))),
        ctx.pow(ctx.num(den), ctx.num(-1)),
    ]
}
```

#### The Method of Undetermined Coefficients

`particular_solution_undetermined` dispatches on the type of the forcing term:

| Forcing type | Detection | Solving strategy |
|---|---|---|
| Polynomial $f(x)$ | `is_polynomial_in` | set $y_p = \sum A_k x^{k+s}$, substitute back, solve for the coefficients |
| Exponential $F e^{kx}$ | `extract_exp_forcing` | set $y_p = A x^s e^{kx}$ |
| Trigonometric $f_c\cos\omega x + f_s\sin\omega x$ | `extract_trig_forcing` | set $y_p = x^s(A\cos\omega x + B\sin\omega x)$ |

$s$ is the resonance shift: when $k$ (or $0$) is an $s$-fold root of the characteristic equation, the trial solution is multiplied by $x^s$.

**Back-substitution algorithm for polynomial forcing**:

For $y_p = \sum_{k=0}^{d} A_k x^{k+s}$, substituting into $ay'' + by' + cy = f(x)$ yields recurrence relations from the coefficients of equal powers of $x$. Solve backwards from the highest power $k = d$:

- $s = 0$ (no resonance): $A_k = \frac{f_k - b(k+1)A_{k+1} - a(k+2)(k+1)A_{k+2}}{c}$
- $s = 1$ (simple resonance): $A_k = \frac{f_k - a(k+2)(k+1)A_{k+1}}{b(k+1)}$
- $s = 2$ (double resonance): $A_k = \frac{f_k}{a(k+2)(k+1)}$

**Superposition principle**: each additive term of the polynomial forcing is solved separately and the results are summed.

#### Variation of Parameters

When the method of undetermined coefficients does not apply, it falls back to variation of parameters (`variation_of_parameters`):

Given the homogeneous solutions $y_1, y_2$ and the standard-form equation $y'' + py' + qy = g$:

1. Compute the Wronskian $W = y_1 y_2' - y_1' y_2$
2. $u_1' = -y_2 g / W$, $u_2' = y_1 g / W$
3. Integrate to obtain $u_1, u_2$
4. $y_p = y_1 u_1 + y_2 u_2$

If the integration returns the unevaluated `Integral(...)` form, the method fails (`is_integral_fallback`).

### Cauchy–Euler Equations

**Form**: $a x^2 y'' + b x y' + c y = f(x)$.

**Algorithm**:

1. Extract the coefficients with `extract_cauchy_euler_coeffs`
2. Reduce to constant coefficients (divide by the powers of $x$)
3. Solve the indicial equation $ar(r-1) + br + c = 0$, i.e. $ar^2 + (b-a)r + c = 0$
4. Construct the homogeneous solution from the discriminant:

| Discriminant | Fundamental solution set |
|---|---|
| $\Delta > 0$ (two real roots $r_1, r_2$) | $x^{r_1},\; x^{r_2}$ |
| $\Delta = 0$ (repeated root $r$) | $x^r,\; x^r \ln x$ |
| $\Delta < 0$ (complex roots $\alpha \pm \beta i$) | $x^\alpha\cos(\beta\ln x),\; x^\alpha\sin(\beta\ln x)$ |

5. For the inhomogeneous term, find a particular solution by variation of parameters (divide by $ax^2$ to get standard form)

### Reduction of Order

**Applies to**: any second-order linear ODE (constant coefficients not required).

**Algorithm**:

1. Try candidate solutions $y_1 \in \{1, x, x^2, e^x, e^{-x}, e^{2x}\}$
2. Verify $a y_1'' + b y_1' + c y_1 = 0$ (`satisfies_extracted`)
3. If $y_1$ satisfies the homogeneous equation, the second solution is:
$$y_2 = y_1 \int \frac{e^{-\int p\,dx}}{y_1^2}\,dx$$
where $p = b/a$
4. The particular solution is obtained by variation of parameters
5. Return $y = C_1 y_1 + C_2 y_2 + y_p$

**Design choice**: the candidate set $\{1, x, x^2, e^x, e^{-x}, e^{2x}\}$ covers the simple solution forms most common in teaching and engineering. For more general equations, the fallback is `PowerSeries`.

---

## The Laplace Transform Method

`laplace.rs` implements Laplace-transform solving of initial-value problems, invoked through the `dsolve_ivp` entry point.

### Applicability

- First order: $a y' + b y = f(x)$, $y(0) = y_0$
- Second order: $a y'' + b y' + c y = f(x)$, $y(0) = y_0$, $y'(0) = y_1$
- The coefficients $a, b, c$ must be integer constants
- The forcing term $f(x)$ must be in the transform table

### The Forward Transform

`laplace_kernel` transforms the $x$-dependent kernel functions. The supported kernels:

| Kernel | Laplace transform |
|---|---|
| $1$ (constant) | $1/s$ |
| $x$ | $1/s^2$ |
| $x^n$ ($n \leq 12$) | $n!/s^{n+1}$ |
| $e^{kx}$ | $1/(s-k)$ |
| $\sin\omega x$ | $\omega/(s^2 + \omega^2)$ |
| $\cos\omega x$ | $s/(s^2 + \omega^2)$ |
| $e^{kx}\sin\omega x$ | $\omega/((s-k)^2 + \omega^2)$ |
| $e^{kx}\cos\omega x$ | $(s-k)/((s-k)^2 + \omega^2)$ |

`split_const_factors` splits a term into an $x$-independent constant and an $x$-dependent kernel, transforming only the kernel part. `linear_coeff_of` extracts $k$ from $kx$ (for $x$ itself it returns $1$).

### Building the Algebraic Equation

The transformed algebraic equation:

**First order**: $a(sY - y_0) + bY = F(s)$
$$Y(s) = \frac{F(s) + a\,y_0}{as + b}$$

**Second order**: $a(s^2Y - sy_0 - y_1) + b(sY - y_0) + cY = F(s)$
$$Y(s) = \frac{F(s) + a(sy_0 + y_1) + by_0}{as^2 + bs + c}$$

### The Inverse Transform

`inverse_laplace_term` performs partial fraction decomposition and table-lookup inversion of the rational function $Y(s) = \frac{n_1 s + n_0}{as^2 + bs + c}$.

First compute the discriminant $\Delta = b^2 - 4ac$; there are three cases depending on the roots:

**Case 1: two distinct real roots $r_1 \neq r_2$**

$$\frac{n_1 s + n_0}{(s - r_1)(s - r_2)} = \frac{A}{s - r_1} + \frac{B}{s - r_2}$$

with $A = \frac{n_1 r_1 + n_0}{r_1 - r_2}$ and $B = \frac{n_1 r_2 + n_0}{r_2 - r_1}$.

Inverse transform: $A\,e^{r_1 x} + B\,e^{r_2 x}$.

**Case 2: repeated root $r$** ($\Delta = 0$)

$$\frac{n_1 s + n_0}{(s - r)^2} = \frac{n_1}{s - r} + \frac{n_1 r + n_0}{(s - r)^2}$$

Inverse transform: $n_1\,e^{rx} + (n_1 r + n_0)\,x\,e^{rx}$.

**Case 3: complex roots $k \pm i\omega$** ($\Delta < 0$)

$$\frac{n_1 s + n_0}{(s - k)^2 + \omega^2} = e^{kx}\!\left[n_1\cos\omega x + \frac{n_0 + kn_1}{\omega}\sin\omega x\right]$$

The $e^{kx}$ factor is dropped when $k = 0$. Requires $(n_0 + kn_1) \bmod \omega = 0$ (integer arithmetic).

### Integer-Arithmetic Constraints

The Laplace module uses `i64` integer arithmetic throughout (via `const_i64`, `quadratic_coeffs_i64`, `linear_coeffs_i64`). This restricts the applicability — equations whose coefficients contain symbolic constants or irrational numbers fall back to other methods. Design rationale: integer arithmetic is completely exact and introduces no floating-point error.

---

## Power Series and the Frobenius Method

`series.rs` implements the two series methods, used as the fallback when every other method fails.

### The Power Series Method

**Applies to**: solutions of linear ODEs near an ordinary point $x_0$.

**Algorithm** (`solve_power_series`):

1. Construct symbolic coefficients $a_0, a_1, \ldots, a_{N-1}$
2. Manually construct the series $S = \sum a_n (x - x_0)^n$ together with its first and second derivatives $S', S''$ (`build_series_triple` — `diff` is not used, to avoid differentiating the $a_n$)
3. Substitute into the ODE: $R(x) = a\,S'' + b\,S' + c\,S - f = 0$
4. Satisfy order by order: for $k = 0, 1, \ldots$, compute $R^{(k)}(x_0) = 0$
5. Each condition determines one new coefficient $a_{k+\text{order}}$ linearly
6. Solved coefficients are substituted back immediately

**Key implementation details**:

- `substitute_series` recursively replaces all $y(x)$, $y'(x)$, $y''(x)$ in the equation by the series forms
- `solve_linear_coeff` solves $a_n$ from the condition equation $\text{coeff} \cdot a_n + \text{rest} = 0$
- The free parameters $a_0, \ldots, a_{\text{order}-1}$ remain symbolic (initial conditions)
- If a condition contains no unknown coefficient but is nonzero, $x_0$ is not an ordinary point → return `None` and hand over to Frobenius

**Limitations**: by default $N = 8$ terms. The coefficients are rational expressions in the initial values (not numerical).

### The Frobenius Method

**Applies to**: solutions of second-order linear ODEs near a regular singular point $x_0 = 0$.

**Algorithm** (`solve_frobenius`):

1. Extract the coefficients $a(x), b(x), c(x)$ and verify $f = 0$ (homogeneous equations only)
2. Set $y = x^r \sum a_n x^n$ and construct $y, y', y''$
3. Substitute into the ODE and group the residual by powers of $x$
4. The lowest-power group yields the **indicial equation** $Ar^2 + Br + C = 0$
5. Find rational roots of the indicial equation (`indicial_coeffs` → `isqrt` + discriminant)
6. Take the larger root $r_1$ and recurse through the subsequent groups to determine the $a_n$

**Construction of $y$**:

$$y = u \cdot S, \quad y' = r x^{-1} u S + u S', \quad y'' = r(r-1)x^{-2}uS + 2rx^{-1}uS' + uS''$$

where $u = x^r$ is an opaque placeholder and $S = \sum a_n x^n$. `strip_x_and_u` strips the powers of $x$ and the $u$ factor from the residual terms, which is used for grouping.

**Extraction of the indicial equation** (`indicial_coeffs`):

The terms of the lowest-power group should be $a_0$ times a quadratic in $r$. Extract $A, B, C$ so that the group equals $a_0(Ar^2 + Br + C)$. `expand_indicial` handles expansions containing the factor $(r-1)$.

**Limitations**:
- Only $x_0 = 0$ is supported
- Homogeneous equations only
- Only real rational roots
- Only the series for the larger root is taken (the smaller root may give a linearly independent solution)

---

## $2\times2$ Linear ODE Systems

`systems.rs` implements solving of constant-coefficient linear systems, invoked through the `dsolve_system` entry point.

### Input Format

The system $\mathbf{Y}' = A\mathbf{Y}$ is given as two equations:
$$\text{Derivative}(y_i, x) - (a_{i1}y_1 + a_{i2}y_2) = 0$$

`extract_coeff` extracts the $2\times2$ coefficient matrix $A$ from each equation.

### The Characteristic Polynomial

$$\lambda^2 - \text{tr}(A)\lambda + \det(A) = 0$$

with $\text{tr} = a_{11} + a_{22}$ and $\det = a_{11}a_{22} - a_{12}a_{21}$.

The discriminant is $\Delta = \text{tr}^2 - 4\det$.

### Solving the Three Cases

#### Distinct Real Eigenvalues

$\lambda_1 \neq \lambda_2$, both integers.

1. `eigenvector` solves $(A - \lambda I)\mathbf{v} = 0$: if $a_{12} \neq 0$, take $\mathbf{v} = (a_{12},\; \lambda - a_{11})$
2. General solution:
$$\mathbf{Y} = C_1 \mathbf{v}_1 e^{\lambda_1 x} + C_2 \mathbf{v}_2 e^{\lambda_2 x}$$

#### Repeated Eigenvalues

$\lambda_1 = \lambda_2 = \lambda$.

**Subcase A**: $A = \lambda I$ (complete matrix):
$$\mathbf{Y} = C_1 \mathbf{e}_1 e^{\lambda x} + C_2 \mathbf{e}_2 e^{\lambda x}$$

**Subcase B**: $A$ is defective:
1. Find the eigenvector $\mathbf{v}$ (`eigenvector`)
2. Find a generalized eigenvector $\mathbf{w}$ satisfying $(A - \lambda I)\mathbf{w} = \mathbf{v}$ (`generalized_eigenvector`)
3. General solution:
$$\mathbf{Y} = e^{\lambda x}\!\left[C_1 \mathbf{v} + C_2(x\mathbf{v} + \mathbf{w})\right]$$

`generalized_eigenvector` solves the linear system $(a_{11}-\lambda)w_1 + a_{12}w_2 = v_1$ for integer solutions, trying $w_1 = 0$ and $w_1 = 1$.

#### Complex Eigenvalues

$\lambda = \alpha \pm \beta i$ (requires $\text{tr}$ even and $-\Delta$ a perfect square).

1. $\alpha = \text{tr}/2$, $\beta = \sqrt{-\Delta}/2$
2. Real part $p$ and imaginary part $q$ of the eigenvector $(v + iw)$: take $v_1 = a_{12}$, $v_2 = (\alpha - a_{11}) + i\beta$, hence $p = (a_{12},\; \alpha - a_{11})$ and $q = (0,\; \beta)$
3. Real fundamental solutions:
$$\mathbf{Y}_1 = e^{\alpha x}(p\cos\beta x - q\sin\beta x)$$
$$\mathbf{Y}_2 = e^{\alpha x}(p\sin\beta x + q\cos\beta x)$$

The $e^{\alpha x}$ factor is dropped when $\alpha = 0$, to avoid a redundant `exp(0)`.

---

## Helper Utilities

`util.rs` provides the infrastructure shared by the ODE module.

### Order Detection

`ode_order` recursively scans the AST for `Derivative(func, var, var, ...)` nodes and returns the highest order of differentiation.

### Linearity Checks

`is_linear_in` checks that `func` and its derivatives occur only to the first power in the expression:
- `is_func_first_degree` counts the occurrences of `func`/`Derivative(func, ...)` inside multiplicative factors
- A factor containing `Pow(func, Num(n))` with $n \geq 2$ means nonlinearity

### Collecting Like Terms

`collect_terms` solves the problem that the simplifier does not collect like terms:

1. Decompose additive terms into a rational coefficient plus base-exponent factor pairs
2. Group by base-exponent signature
3. Sum the rational coefficients within each group
4. Drop groups whose coefficient is zero

The decomposition supports rational exponents: $x \cdot x^{-1/2} = x^{1/2}$.

### Exponential Simplification

`exp_simplify` handles $\exp(k \cdot \log u) = u^k$:
- Recognizes exponents of the form `Mul(k, Fun("log", [u]))`
- Numerical constant factors are discarded from the $\log$ argument (no effect for ODE integrating factors)
- Falls back to the literal `exp` form

### Substituting and Verifying Solutions

`substitute_solution` substitutes a candidate solution $y = \text{sol}$ into the equation:
- Replace `func` → `sol`
- Replace `Derivative(func, var)` → `diff(sol, var)`
- Replace `Derivative(func, var, var)` → `diff(diff(sol, var), var)`

`verify_solution` (tests only) substitutes and checks whether the residual simplifies to zero.

---

## Design Choices and Trade-offs

### Simplification Strategy

The ODE module relies heavily on `ocas_rewrite::simplify` (rule-based simplification, at most 20 rounds of iteration) and `ocas_atom::normalize` (structural normalization). Their division of labor:

- `normalize`: fast structural equivalence testing (string comparison)
- `simplify`: rule-based algebraic simplification (combining, expanding, eliminating)

In exactness detection and solution construction, the two are used alternately to handle equivalences the simplifier cannot discover automatically.

### Integration Dependency

The ODE solvers call `crate::integral::integrate` to compute indefinite integrals. This launches the full integration pipeline (fast table → rational function → Risch → trigonometric rewrite → special functions → heuristic techniques → unevaluated form). When the integral cannot be evaluated in closed form, the `Integral(...)` form is returned — in that case variation of parameters and reduction of order are marked as failed.

### Integer Restrictions

The Laplace transform module strictly uses `i64` integer arithmetic and rejects symbolic coefficients. This is because partial fraction decomposition requires exact rational arithmetic, and integers are the safest choice over `i64`. For initial-value problems requiring symbolic coefficients, one should fall back to the general `dsolve` plus substitution of the initial conditions.

### The Fallback Chain

The priority order returned by the classifier ensures the most specific method is tried first:

```
LinearFirst → Bernoulli → Separable → Exact → Homogeneous
    → LinearConstantCoeff → CauchyEuler → ReductionOfOrder
        → PowerSeries
```

`PowerSeries` is the method of last resort and always appears at the end of
the list; its solver tries the power-series method first and falls back to
Frobenius internally. `dsolve` tries them in order until the first success.

---

## Source Locations

| File | Responsibility |
|---|---|
| `ocas-calc/src/ode/mod.rs` | Entry points: `dsolve`, `dsolve_ivp`, `dsolve_system`, `ODE`, `ODESolution` |
| `ocas-calc/src/ode/classify.rs` | Classifier: `classify_ode`, `ODEType`, the detection functions |
| `ocas-calc/src/ode/first_order.rs` | First-order solvers: separable, linear, Bernoulli, exact, homogeneous |
| `ocas-calc/src/ode/second_order.rs` | Second-order solvers: constant coefficient, Cauchy–Euler, reduction of order, variation of parameters |
| `ocas-calc/src/ode/laplace.rs` | Laplace transform: forward transform table, inverse transform (three root types), IVP solving |
| `ocas-calc/src/ode/series.rs` | Power series and Frobenius methods |
| `ocas-calc/src/ode/systems.rs` | $2\times2$ systems: eigenvalue decomposition, real/complex/repeated-root cases |
| `ocas-calc/src/ode/util.rs` | Utilities: order, linearity, like-term collection, exponential simplification, solution substitution |

---

## See Also

- [ODE Solving Theory](../math/ode-theory.md) — mathematical derivations and proofs of each method
- [Symbolic Calculus](../math/symbolic-calculus.md) — theoretical foundations of differentiation and Taylor expansion
- [Rust API: Solvers](../api/rust-solvers.md) — API reference for `dsolve`, `dsolve_ivp`, `dsolve_system`
- [Rust API: Calculus](../api/rust-calculus.md) — API reference for `diff`, `integrate`
