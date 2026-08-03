# Advanced: ODE Solving Theory

## Prerequisites

- [Polynomial Algebra](./polynomial-algebra.md) — polynomial operations and factorization
- [Linear Algebra](./linear-algebra.md) — matrix operations and eigenvalues
- [Symbolic Calculus](./symbolic-calculus.md) — the basic rules of differentiation and integration

---

## Basic Concepts

### Definition of an Ordinary Differential Equation

An **ordinary differential equation** (ODE) is an equation involving an unknown function $y(x)$ and its derivatives. The general form is

$$F\!\left(x,\, y,\, y',\, y'',\, \dots,\, y^{(n)}\right) = 0$$

where $y^{(k)} = \dfrac{d^k y}{dx^k}$.

### Order

The **highest-order derivative** appearing in the ODE determines the order of the equation. For example:

| Equation | Order |
|---|---|
| $y' + y = 0$ | 1 |
| $y'' + 2y' + y = \sin x$ | 2 |
| $y''' - y = 0$ | 3 |

### Linear and Nonlinear

An $n$-th order ODE is **linear** if and only if it can be written as

$$a_n(x)\, y^{(n)} + a_{n-1}(x)\, y^{(n-1)} + \cdots + a_1(x)\, y' + a_0(x)\, y = g(x)$$

where the $a_k(x)$ and $g(x)$ depend only on the independent variable $x$, not on $y$. If $g(x) \equiv 0$, the equation is called **homogeneous** linear; otherwise it is **non-homogeneous**.

Equations containing terms like $y^2$, $\sin y$, or $y\,y'$ are **nonlinear**.

### Initial Conditions and Boundary Conditions

- **Initial value problem** (IVP): given $y(x_0) = y_0$ (first order), or $y(x_0) = y_0,\; y'(x_0) = y_1$ (second order), the solution is uniquely determined.
- **Boundary value problem** (BVP): conditions are specified at two different points, e.g. $y(0) = 0,\; y(\pi) = 0$.

### The Existence and Uniqueness Theorem (Picard–Lindelöf)

**Theorem** (first-order case). Let $f(x, y)$ be continuous on the rectangle $R = \{(x,y) : |x - x_0| \leq a,\; |y - y_0| \leq b\}$ and satisfy the Lipschitz condition in $y$:

$$|f(x, y_1) - f(x, y_2)| \leq L\, |y_1 - y_2|$$

Then the initial value problem $y' = f(x,y),\; y(x_0) = y_0$ has a **unique** solution on $|x - x_0| \leq h$ ($h = \min(a,\, b/M)$, $M = \max_R |f|$).

**Generalization to order $n$**. Convert the $n$-th order equation into a first-order system $\mathbf{Y}' = \mathbf{F}(x, \mathbf{Y})$ and impose the same Lipschitz condition on $\mathbf{F}$.

This theorem guarantees that when the right-hand side is sufficiently "nice", the solution of an ODE exists and is unique. It is the theoretical foundation of all symbolic solvers — if a symbolic method finds a solution, Picard–Lindelöf guarantees it is *the* unique solution.

---

## Core Theory

### First-Order ODEs

#### Separable Equations

**Form**:

$$g(y)\, y' = f(x)$$

**Solution steps**:

1. Write the equation as $g(y)\, dy = f(x)\, dx$
2. Integrate both sides: $\displaystyle\int g(y)\, dy = \int f(x)\, dx + C$

**Example**. Solve $y' = xy$:

$$\frac{dy}{y} = x\, dx \quad\Longrightarrow\quad \ln|y| = \frac{x^2}{2} + C \quad\Longrightarrow\quad y = C_1\, e^{x^2/2}$$

**Implementation in oCAS**. The classifier `is_separable` decides separability by checking whether the equation can be written so that "all terms containing $y$ are paired with $y'$, and the remaining terms contain only $x$". The solver `solve_separable` splits the equation into a $(y\text{-part},\, x\text{-part})$ pair, integrates each part, and constructs the implicit solution $\int g\, dy - \int f\, dx = C$.

---

#### Linear First-Order Equations

**Form**:

$$y' + p(x)\, y = q(x)$$

**Integrating factor method**. Multiply by the integrating factor

$$\mu(x) = e^{\int p(x)\, dx}$$

The equation becomes

$$\frac{d}{dx}\!\bigl[\mu(x)\, y\bigr] = \mu(x)\, q(x)$$

Integrating gives

$$y = \frac{1}{\mu(x)} \int \mu(x)\, q(x)\, dx + \frac{C}{\mu(x)}$$

**Derivation**. Observe that $\mu' = p\,\mu$, hence

$$\mu\, y' + p\,\mu\, y = \mu\, y' + \mu'\, y = (\mu\, y)'$$

**Example**. Solve $y' - \dfrac{2}{x}\, y = x^2$:

$$p(x) = -\frac{2}{x}, \quad \mu = e^{\int -2/x\, dx} = e^{-2\ln x} = x^{-2}$$

$$\frac{d}{dx}\!\bigl[x^{-2} y\bigr] = x^{-2} \cdot x^2 = 1$$

$$x^{-2}\, y = x + C \quad\Longrightarrow\quad y = x^3 + Cx^2$$

**Implementation in oCAS**. `solve_linear_first` extracts $p(x)$ and $q(x)$, computes $\mu = \exp\!\bigl(\int p\, dx\bigr)$, and returns

$$y = \mu^{-1}\!\left(\int \mu\, q\, dx + C\right)$$

---

#### Bernoulli Equations

**Form**:

$$y' + p(x)\, y = q(x)\, y^n, \qquad n \neq 0,\, 1$$

**Transformation**. Let $v = y^{1-n}$; then $v' = (1-n)\, y^{-n}\, y'$. Multiplying both sides of the original equation by $(1-n)\, y^{-n}$:

$$(1-n)\, y^{-n}\, y' + (1-n)\, p(x)\, y^{1-n} = (1-n)\, q(x)$$

i.e.

$$v' + (1-n)\, p(x)\, v = (1-n)\, q(x)$$

This is a linear first-order equation in $v$; solve it with the integrating factor method and substitute back $v = y^{1-n}$.

**Example**. Solve $y' + \dfrac{1}{x}\, y = x\, y^2$ ($n = 2$):

Let $v = y^{-1}$, so $v' = -y^{-2}\, y'$, giving

$$-v' + \frac{1}{x}\, v = x \quad\Longrightarrow\quad v' - \frac{1}{x}\, v = -x$$

With the integrating factor $\mu = x^{-1}$, we get $v = -x^2 + Cx$, i.e. $y = \dfrac{1}{Cx - x^2}$.

**Implementation in oCAS**. `solve_bernoulli` detects the nonlinear power $n$, constructs the substitution $v = y^{1-n}$, calls `solve_linear_first` on the linearized equation, and then substitutes back.

---

#### Exact Equations

**Form**:

$$M(x,y) + N(x,y)\, y' = 0$$

If there exists a potential function $F(x,y)$ with $\dfrac{\partial F}{\partial x} = M$ and $\dfrac{\partial F}{\partial y} = N$, then the equation is equivalent to $\dfrac{d}{dx} F(x, y(x)) = 0$, and the solution is $F(x,y) = C$.

**Exactness criterion**. The equation is exact if and only if

$$\frac{\partial M}{\partial y} = \frac{\partial N}{\partial x}$$

**Solution steps**:

1. Integrate $\partial F/\partial x = M$ with respect to $x$: $F = \int M\, dx + h(y)$
2. Determine $h(y)$ from $\partial F/\partial y = N$: $h'(y) = N - \partial/\partial y \int M\, dx$
3. Substitute to obtain $F(x,y) = C$

**Integrating factors for non-exact equations**. If the equation is not exact, one may try to find an integrating factor $\mu$:

- if $\dfrac{M_y - N_x}{N}$ depends only on $x$, then $\mu(x) = \exp\!\left(\int \dfrac{M_y - N_x}{N}\, dx\right)$
- if $\dfrac{N_x - M_y}{M}$ depends only on $y$, then $\mu(y) = \exp\!\left(\int \dfrac{N_x - M_y}{M}\, dy\right)$

**Example**. Solve $(2xy + 3) + (x^2 + 4y)\, y' = 0$:

$M_y = 2x = N_x$, so the equation is exact.

$$F = \int (2xy + 3)\, dx = x^2 y + 3x + h(y)$$

$$F_y = x^2 + h'(y) = x^2 + 4y \quad\Longrightarrow\quad h'(y) = 4y \quad\Longrightarrow\quad h = 2y^2$$

Solution: $x^2 y + 3x + 2y^2 = C$.

**Implementation in oCAS**. `solve_exact` first verifies exactness with `partials_equal`. If the equation is not exact, `find_integrating_factor` tries the two integrating factors above (the cases depending only on $x$ or only on $y$). Once the potential function is found, an implicit solution is returned.

---

#### Homogeneous Equations

**Form**:

$$y' = f\!\left(\frac{y}{x}\right)$$

Equivalently, the right-hand side $f(x,y)$ is a degree-zero homogeneous function of $(x,y)$: $f(tx, ty) = f(x,y)$.

**Transformation**. Let $v = y/x$, so $y = vx$ and $y' = v + xv'$. Substituting:

$$v + x\, v' = f(v) \quad\Longrightarrow\quad x\, v' = f(v) - v$$

This is a separable equation in $v$ and $x$:

$$\frac{dv}{f(v) - v} = \frac{dx}{x}$$

Integrate and substitute back $v = y/x$.

**Example**. Solve $y' = \dfrac{y + x}{x} = \dfrac{y}{x} + 1$:

Let $v = y/x$; then $v + xv' = v + 1$, so $xv' = 1$, giving $v = \ln|x| + C$ and $y = x\ln|x| + Cx$.

**Implementation in oCAS**. The classifier `is_homogeneous` requires the equation to be first-order **linear**, with all additive terms having the same total degree in $(x, y)$ (degree-zero homogeneity), and with a term free of $y$ present. The solver `solve_homogeneous` performs the substitution $v = y/x$, reducing the equation to a separable one.

---

### Second-Order Linear ODEs

The general form of a second-order linear ODE is

$$a(x)\, y'' + b(x)\, y' + c(x)\, y = f(x)$$

**Superposition principle**. If $y_1, y_2$ are two linearly independent solutions of the homogeneous equation and $y_p$ is a particular solution of the non-homogeneous equation, then the general solution is

$$y = C_1\, y_1 + C_2\, y_2 + y_p$$

#### Constant Coefficient Equations

**Form**:

$$a\, y'' + b\, y' + c\, y = f(x), \qquad a, b, c \in \mathbb{R},\; a \neq 0$$

**Characteristic equation**. The characteristic equation of the homogeneous equation $ay'' + by' + cy = 0$ is

$$a\, r^2 + b\, r + c = 0$$

The discriminant $\Delta = b^2 - 4ac$ determines the form of the solution:

| Discriminant | Characteristic roots | General homogeneous solution |
|---|---|---|
| $\Delta > 0$ | $r_1, r_2 \in \mathbb{R},\; r_1 \neq r_2$ | $C_1\, e^{r_1 x} + C_2\, e^{r_2 x}$ |
| $\Delta = 0$ | $r_1 = r_2 = r$ | $(C_1 + C_2\, x)\, e^{rx}$ |
| $\Delta < 0$ | $\alpha \pm \beta i$ | $e^{\alpha x}(C_1 \cos\beta x + C_2 \sin\beta x)$ |

**The method of undetermined coefficients**. When $f(x)$ is a polynomial, an exponential $e^{kx}$, a trigonometric function $\sin\omega x$, $\cos\omega x$, or a product/sum of these, the form of the particular solution can be guessed.

| Form of $f(x)$ | Guessed form of $y_p$ |
|---|---|
| $P_n(x)$ (a polynomial of degree $n$) | $Q_n(x) = a_0 + a_1 x + \cdots + a_n x^n$ |
| $F\, e^{kx}$ | $A\, e^{kx}$ |
| $f_c \cos\omega x + f_s \sin\omega x$ | $A\cos\omega x + B\sin\omega x$ |
| $F\, e^{kx} \cos\omega x$ | $e^{kx}(A\cos\omega x + B\sin\omega x)$ |

**Resonance**. When the guessed form is linearly dependent with the homogeneous solutions, multiply by $x^s$ ($s$ the multiplicity of $k$ as a characteristic root):

- if $k$ is a simple characteristic root: $y_p = Ax\, e^{kx}$
- if $k$ is a double characteristic root: $y_p = Ax^2\, e^{kx}$

For polynomial forcing, if $r = 0$ is a characteristic root of multiplicity $s$, the guessed form must be multiplied by $x^s$.

**Example**. Solve $y'' - 3y' + 2y = e^{3x}$:

The characteristic equation is $r^2 - 3r + 2 = 0$, with roots $r_1 = 1, r_2 = 2$.

General homogeneous solution: $C_1 e^x + C_2 e^{2x}$.

$k = 3$ is not a characteristic root, so guess $y_p = Ae^{3x}$. Substituting: $9A - 9A + 2A = 1$, so $A = 1/2$.

General solution: $y = C_1 e^x + C_2 e^{2x} + \dfrac{1}{2} e^{3x}$.

**Implementation in oCAS**. `solve_constant_coeff` extracts the coefficients $(a, b, c)$ and the forcing term $f(x)$. `constant_coeff_basis` constructs the homogeneous basis solutions according to the discriminant. `particular_solution_undetermined` handles polynomial, exponential, and trigonometric forcing, including complete resonance detection.

---

#### Cauchy–Euler Equations

**Form**:

$$a\, x^2\, y'' + b\, x\, y' + c\, y = f(x)$$

The coefficient of $y^{(k)}$ is a constant multiple of $x^k$.

**Transformation**. Let $x = e^t$ (i.e. $t = \ln x$), using

$$y' = \frac{1}{x}\, \frac{dy}{dt}, \qquad y'' = \frac{1}{x^2}\!\left(\frac{d^2 y}{dt^2} - \frac{dy}{dt}\right)$$

The original equation becomes a constant coefficient equation

$$a\, \ddot{y} + (b - a)\, \dot{y} + c\, y = f(e^t)$$

where $\dot{y} = dy/dt$.

**Indicial equation**. For the homogeneous equation, guess $y = x^r$ directly and substitute to obtain the indicial equation

$$a\, r(r-1) + b\, r + c = 0 \quad\Longleftrightarrow\quad a\, r^2 + (b - a)\, r + c = 0$$

The structure of the solutions is analogous to the constant coefficient case (three cases according to the discriminant), but the fundamental solutions are $x^r$ instead of $e^{rx}$:

| Discriminant | Fundamental solutions |
|---|---|
| $\Delta > 0$ | $x^{r_1},\; x^{r_2}$ |
| $\Delta = 0$ | $x^r,\; x^r \ln x$ |
| $\Delta < 0$ ($\alpha \pm \beta i$) | $x^\alpha \cos(\beta\ln x),\; x^\alpha \sin(\beta\ln x)$ |

**Example**. Solve $x^2 y'' - 2xy' + 2y = 0$:

Indicial equation: $r(r-1) - 2r + 2 = r^2 - 3r + 2 = 0$, with roots $r_1 = 1, r_2 = 2$.

General solution: $y = C_1 x + C_2 x^2$.

**Implementation in oCAS**. The classifier `is_cauchy_euler` checks whether the equation matches the pattern $c_k x^k y^{(k)}$. `solve_cauchy_euler` extracts the Cauchy–Euler coefficients, solves the indicial equation, and constructs the general solution according to the discriminant. `cauchy_euler_basis` handles the three discriminant cases.

---

#### Reduction of Order

When one solution $y_1$ of the homogeneous equation is known, a second linearly independent solution can be obtained by reduction of order.

**Method**. Let $y_2 = v(x)\, y_1$ and substitute into $y'' + p\, y' + q\, y = 0$. Using the fact that $y_1$ satisfies the homogeneous equation, the equation simplifies to a first-order equation in $v'$:

$$v'' y_1 + v'(2y_1' + p\, y_1) = 0$$

Let $w = v'$ and separate variables:

$$\frac{w'}{w} = -\frac{2y_1' + p\, y_1}{y_1} = -2\, \frac{y_1'}{y_1} - p$$

Integrating gives

$$w = \frac{e^{-\int p\, dx}}{y_1^2}$$

Hence

$$y_2 = y_1 \int \frac{e^{-\int p\, dx}}{y_1^2}\, dx$$

**Example**. $x^2 y'' - 2xy' + 2y = 0$, with known solution $y_1 = x$.

Standard form: $y'' - \dfrac{2}{x}\, y' + \dfrac{2}{x^2}\, y = 0$, so $p = -2/x$.

$$y_2 = x \int \frac{e^{\int 2/x\, dx}}{x^2}\, dx = x \int \frac{x^2}{x^2}\, dx = x \cdot x = x^2$$

**Implementation in oCAS**. `solve_reduction_of_order` first tries simple candidate solutions ($1, x, x^2, e^x, e^{-x}, e^{2x}$) in the homogeneous equation. Once $y_1$ is found, the formula above computes $y_2$. For non-homogeneous equations, variation of parameters is then used to find a particular solution.

---

#### Variation of Parameters

When the method of undetermined coefficients does not apply (e.g. $f(x)$ is not of standard form), variation of parameters is the general method.

**Method**. Given two linearly independent solutions $y_1, y_2$ of the homogeneous equation, look for a particular solution of the non-homogeneous equation of the form

$$y_p = u_1(x)\, y_1 + u_2(x)\, y_2$$

where $u_1, u_2$ satisfy

$$\begin{cases} u_1'\, y_1 + u_2'\, y_2 = 0 \\ u_1'\, y_1' + u_2'\, y_2' = g(x) \end{cases}$$

Here $g(x) = f(x)/a$ is the right-hand side of the standard form. By Cramer's rule:

$$u_1' = -\frac{y_2\, g}{W}, \qquad u_2' = \frac{y_1\, g}{W}$$

where $W = y_1\, y_2' - y_1'\, y_2$ is the **Wronskian**.

Therefore

$$y_p = -y_1 \int \frac{y_2\, g}{W}\, dx + y_2 \int \frac{y_1\, g}{W}\, dx$$

**Significance of the Wronskian**. $y_1, y_2$ are linearly independent if and only if $W \neq 0$. For a second-order homogeneous linear ODE, Abel's formula gives $W(x) = W(x_0)\, \exp\!\bigl(-\int_{x_0}^x p\, dt\bigr)$, so $W$ is either identically zero or never zero.

**Implementation in oCAS**. `variation_of_parameters` takes the two basis solutions $y_1, y_2$ and the standard-form right-hand side $g$, computes the Wronskian, and then uses oCAS's symbolic integrator (`integrate`) to compute $u_1, u_2$. If the integrals cannot be evaluated in closed form (returning an unevaluated `Integral(...)`), the method is abandoned.

---

### Series Solutions

When the coefficients of an ODE are not constant and no elementary method applies, series methods provide a systematic way to construct solutions.

#### Power Series Solutions (Ordinary Points)

**Applicability**. Consider the second-order linear ODE

$$y'' + p(x)\, y' + q(x)\, y = 0$$

If $p(x)$ and $q(x)$ are analytic at $x = x_0$ (i.e. they can be expanded in convergent power series), then $x_0$ is an **ordinary point**, and the equation has two linearly independent analytic solutions near $x_0$.

**Method**. Let

$$y = \sum_{n=0}^{\infty} a_n (x - x_0)^n$$

Then

$$y' = \sum_{n=1}^{\infty} n\, a_n (x - x_0)^{n-1}, \qquad y'' = \sum_{n=2}^{\infty} n(n-1)\, a_n (x - x_0)^{n-2}$$

Substitute the power series of $y, y', y''$ and of $p(x), q(x)$ into the equation, collect like powers, and set the coefficient of each power to zero, obtaining a **recurrence relation** for the $a_n$.

For a $k$-th order ODE, the initial coefficients $a_0, a_1, \dots, a_{k-1}$ are free parameters (corresponding to the $k$ arbitrary constants of the general solution); the remaining coefficients are determined by the recurrence.

**Example**. Solve $y'' - xy = 0$ (the Airy equation) by a series about $x = 0$:

$$\sum_{n=2}^{\infty} n(n-1)a_n x^{n-2} - x \sum_{n=0}^{\infty} a_n x^n = 0$$

$$\sum_{n=0}^{\infty} (n+2)(n+1)a_{n+2} x^n - \sum_{n=0}^{\infty} a_n x^{n+1} = 0$$

The $x^0$ term: $2 \cdot 1 \cdot a_2 = 0$, so $a_2 = 0$.

The $x^n$ ($n \geq 1$) term: $(n+2)(n+1)\, a_{n+2} - a_{n-1} = 0$, i.e.

$$a_{n+2} = \frac{a_{n-1}}{(n+2)(n+1)}$$

Taking $a_0 = 1, a_1 = 0$: $a_3 = a_0/(3 \cdot 2) = 1/6$, $a_6 = a_3/(6 \cdot 5) = 1/180$, …

Taking $a_0 = 0, a_1 = 1$: $a_4 = a_1/(4 \cdot 3) = 1/12$, $a_7 = a_4/(7 \cdot 6) = 1/504$, …

**Implementation in oCAS**. `solve_power_series` sets $y = \sum a_n x^n$, constructs the series triple of $y, y', y''$ (`build_series_triple`), substitutes it into the equation, groups by powers of $x$, and solves the linear equations for the $a_n$ (`solve_linear_coeff`). It returns the truncated series as `ODESolution::Series(expr, n\_terms)`.

---

#### The Frobenius Method (Regular Singular Points)

**Applicability**. If $x = x_0$ is not an ordinary point, but $(x - x_0)\, p(x)$ and $(x - x_0)^2\, q(x)$ are analytic at $x_0$, then $x_0$ is a **regular singular point**.

**Method**. Let

$$y = x^r \sum_{n=0}^{\infty} a_n x^n = \sum_{n=0}^{\infty} a_n x^{n+r}, \qquad a_0 \neq 0$$

Substitute into the equation. The coefficient of the lowest power $x^{r-2}$ (from the Cauchy–Euler-type leading part) gives the **indicial equation**:

$$A\, r(r-1) + B\, r + C = 0$$

where $A, B, C$ are determined by the coefficients of the equation at $x_0$. The higher powers give the recurrence relation for the $a_n$.

**Relations between the indicial roots**. Let the two roots of the indicial equation be $r_1, r_2$ ($r_1 \geq r_2$):

| Roots | Form of the solutions |
|---|---|
| $r_1 - r_2 \notin \mathbb{Z}$ | two independent series solutions $y_1 = x^{r_1}\sum a_n x^n$, $y_2 = x^{r_2}\sum b_n x^n$ |
| $r_1 = r_2$ | $y_1 = x^{r_1}\sum a_n x^n$, $y_2 = y_1 \ln x + x^{r_1}\sum c_n x^n$ |
| $r_1 - r_2 \in \mathbb{Z}^+$ | $y_1$ as above; $y_2$ may contain a $\ln x$ term (the logarithmic case) |

**Example**. Solve $x^2 y'' + xy' + (x^2 - 1)y = 0$ (Bessel's equation, $\nu = 1$) near $x = 0$:

$p(x) = 1/x$ and $q(x) = 1 - 1/x^2$. So $x = 0$ is a regular singular point.

Let $y = \sum a_n x^{n+r}$. After substitution, the coefficient of $x^r$ is

$$a_0[r(r-1) + r - 1] = a_0(r^2 - 1) = 0$$

Indicial equation: $r^2 - 1 = 0$, with roots $r_1 = 1, r_2 = -1$.

Since $r_1 - r_2 = 2 \in \mathbb{Z}$, the logarithmic case must be checked.

**Implementation in oCAS**. `solve_frobenius` sets $y = x^r \sum a_n x^n$ ($r$ a symbolic variable), substitutes into the equation, and groups by powers of $x$. The lowest-power group is used by `indicial_coeffs` to extract the integer coefficients of the indicial equation $Ar^2 + Br + C = 0$. After solving for $r$, the higher groups determine the $a_n$ by recurrence. Currently only the case of **rational real roots** of the indicial equation is handled, returning the series corresponding to the larger root.

---

### The Laplace Transform

The Laplace transform is a powerful tool for solving initial value problems for linear constant-coefficient ODEs, especially with piecewise continuous forcing or impulse functions.

#### Definition

The Laplace transform of a function $f(x)$ ($x \geq 0$) is defined as

$$\mathcal{L}\{f\}(s) = F(s) = \int_0^{\infty} e^{-sx}\, f(x)\, dx$$

#### Transforms of Derivatives

By integration by parts:

$$\mathcal{L}\{y'\}(s) = s\, Y(s) - y(0)$$

$$\mathcal{L}\{y''\}(s) = s^2\, Y(s) - s\, y(0) - y'(0)$$

In general

$$\mathcal{L}\{y^{(n)}\}(s) = s^n\, Y(s) - \sum_{k=0}^{n-1} s^{n-1-k}\, y^{(k)}(0)$$

#### Common Transform Pairs

| $f(x)$ | $F(s) = \mathcal{L}\{f\}(s)$ |
|---|---|
| $1$ | $\dfrac{1}{s}$ |
| $x^n$ | $\dfrac{n!}{s^{n+1}}$ |
| $e^{kx}$ | $\dfrac{1}{s - k}$ |
| $\cos\omega x$ | $\dfrac{s}{s^2 + \omega^2}$ |
| $\sin\omega x$ | $\dfrac{\omega}{s^2 + \omega^2}$ |
| $e^{kx}\cos\omega x$ | $\dfrac{s - k}{(s-k)^2 + \omega^2}$ |
| $e^{kx}\sin\omega x$ | $\dfrac{\omega}{(s-k)^2 + \omega^2}$ |
| $x\, e^{kx}$ | $\dfrac{1}{(s-k)^2}$ |

#### Steps for Solving an IVP

Take the second-order equation $ay'' + by' + cy = f(x)$, $y(0) = y_0$, $y'(0) = y_1$ as an example:

1. **Transform**: apply the Laplace transform to both sides of the equation

$$a[s^2 Y - sy_0 - y_1] + b[sY - y_0] + cY = F(s)$$

2. **Solve the algebraic equation**:

$$Y(s) = \frac{F(s) + a(sy_0 + y_1) + by_0}{as^2 + bs + c}$$

3. **Partial fraction decomposition**: decompose $Y(s)$ into a linear combination of standard transform pairs

4. **Inverse transform**: look up the table to obtain $y(x) = \mathcal{L}^{-1}\{Y(s)\}(x)$

**Key techniques for the inverse transform**. For rational functions with a quadratic denominator:

- **Distinct real roots** $r_1 \neq r_2$: $\dfrac{n_1 s + n_0}{(s-r_1)(s-r_2)} = \dfrac{A}{s-r_1} + \dfrac{B}{s-r_2}$, with $A = \dfrac{n_1 r_1 + n_0}{r_1 - r_2}$

- **Repeated root** $r$: $\dfrac{n_1 s + n_0}{(s-r)^2} = \dfrac{n_1}{s-r} + \dfrac{n_1 r + n_0}{(s-r)^2}$

- **Complex conjugate roots** $k \pm i\omega$: $\dfrac{n_1 s + n_0}{(s-k)^2 + \omega^2} = e^{kx}\!\left[n_1 \cos\omega x + \dfrac{n_0 + kn_1}{\omega}\sin\omega x\right]$

**Example**. Solve $y'' - 3y' + 2y = 4$, $y(0) = 2$, $y'(0) = 3$:

$$s^2 Y - 2s - 3 - 3(sY - 2) + 2Y = \frac{4}{s}$$

$$(s^2 - 3s + 2)Y = \frac{4}{s} + 2s - 3 = \frac{2s^2 - 3s + 4}{s}$$

$$Y = \frac{2s^2 - 3s + 4}{s(s-1)(s-2)}$$

Partial fractions: $Y = \dfrac{2}{s} - \dfrac{3}{s-1} + \dfrac{3}{s-2}$.

Inverse transform: $y(x) = 2 - 3e^x + 3e^{2x}$.

Verification: $y(0) = 2$, $y'(0) = -3 + 6 = 3$, and $y'' - 3y' + 2y = 4$ ✓.

**Implementation in oCAS**. `dsolve_ivp` is the entry point for IVP solving; internally it calls `solve_laplace`.

The workflow of `solve_laplace`:

1. Determine whether the equation is a first- or second-order linear constant-coefficient equation
2. Extract the coefficients $(a, b, c)$ and the forcing term $f(x)$
3. `laplace_transform` computes $\mathcal{L}\{f\}(s)$ — supporting polynomials, exponentials, trigonometric functions, and their products
4. Construct $Y(s) = \dfrac{F(s) + \text{initial value terms}}{as^2 + bs + c}$
5. `inverse_laplace` decomposes the rational function with `split_fraction` and calls `inverse_distinct_roots`, `inverse_repeated_root`, or `inverse_complex_roots` on each term
6. Return `ODESolution::Explicit(y(x))` (without arbitrary constants)

---

### 2×2 Linear Systems with Constant Coefficients

Consider a first-order linear system with constant coefficients

$$\mathbf{Y}' = A\, \mathbf{Y}, \qquad \mathbf{Y} = \begin{pmatrix} y_1 \\ y_2 \end{pmatrix}, \quad A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}$$

#### General Solution: the Matrix Exponential

The formal solution is $\mathbf{Y}(x) = e^{Ax}\, \mathbf{Y}(0)$. Practical solving depends on the eigenvalues of $A$.

#### Eigenvalue Decomposition

The characteristic equation of $A$ is

$$\lambda^2 - (a_{11} + a_{22})\lambda + (a_{11}a_{22} - a_{12}a_{21}) = 0$$

Let the eigenvalues be $\lambda_1, \lambda_2$ with corresponding eigenvectors $\mathbf{v}_1, \mathbf{v}_2$.

**Case 1: distinct real eigenvalues** ($\lambda_1 \neq \lambda_2 \in \mathbb{R}$)

$$\mathbf{Y} = C_1\, e^{\lambda_1 x}\, \mathbf{v}_1 + C_2\, e^{\lambda_2 x}\, \mathbf{v}_2$$

**Case 2: repeated eigenvalue** ($\lambda_1 = \lambda_2 = \lambda$)

If $A$ is not diagonalizable (only one linearly independent eigenvector $\mathbf{v}$), a generalized eigenvector $\mathbf{w}$ satisfying $(A - \lambda I)\mathbf{w} = \mathbf{v}$ is needed:

$$\mathbf{Y} = C_1\, e^{\lambda x}\, \mathbf{v} + C_2\, e^{\lambda x}(x\, \mathbf{v} + \mathbf{w})$$

**Case 3: complex conjugate eigenvalues** ($\lambda = \alpha \pm \beta i$, $\beta \neq 0$)

The eigenvector also has real and imaginary parts: $\mathbf{v} = \mathbf{p} + i\mathbf{q}$. The real-valued general solution is

$$\mathbf{Y} = e^{\alpha x}\bigl[C_1(\mathbf{p}\cos\beta x - \mathbf{q}\sin\beta x) + C_2(\mathbf{p}\sin\beta x + \mathbf{q}\cos\beta x)\bigr]$$

**Example**. Solve the harmonic oscillator system $y_1' = y_2,\; y_2' = -y_1$:

$$A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$$

Characteristic equation: $\lambda^2 + 1 = 0$, so $\lambda = \pm i$ ($\alpha = 0, \beta = 1$).

For $\lambda = i$: $(A - iI)\mathbf{v} = 0$, giving $\mathbf{v} = (1, i)^T$, $\mathbf{p} = (1, 0)^T$, $\mathbf{q} = (0, 1)^T$.

$$\mathbf{Y} = C_1 \begin{pmatrix} \cos x \\ -\sin x \end{pmatrix} + C_2 \begin{pmatrix} \sin x \\ \cos x \end{pmatrix}$$

i.e. $y_1 = C_1\cos x + C_2\sin x$, $y_2 = -C_1\sin x + C_2\cos x$.

**Implementation in oCAS**. `solve_linear_system` handles $2 \times 2$ constant-coefficient systems:

1. Extract the entries of the matrix $A$ from the system (`extract_coeff`)
2. Compute the characteristic polynomial and find the eigenvalues (within the integers or rationals)
3. For distinct real eigenvalues: `eigenvector` finds integer eigenvectors
4. For a repeated eigenvalue: `generalized_eigenvector` finds the generalized eigenvector $(A - \lambda I)\mathbf{w} = \mathbf{v}$
5. For complex eigenvalues: `solve_complex_2x2` constructs the real-valued fundamental solutions
6. Return `ODESolution::System([y_1, y_2])`, containing the arbitrary constants $C_1, C_2$

---

## Implementation in oCAS

### The Classifier (classify.rs)

`classify_ode` is the entry-point decision maker of ODE solving. It analyzes the structure of the equation and returns the list of candidate methods in priority order:

| Priority | `ODEType` | Detection condition |
|---|---|---|
| 1 | `LinearFirst` | first order; $y$ and $y'$ appear linearly, no nonlinear terms like $y^2, \sin y$, and $y$ is present |
| 2 | `Bernoulli` | first order; contains a nonlinear power $y^n$ ($n \geq 2$) |
| 3 | `Separable` | first order; the terms containing $y$ can be grouped with $y'$, the remaining terms contain only $x$ |
| 4 | `Exact` | first order; $M + Ny' = 0$ with $\partial M/\partial y = \partial N/\partial x$ |
| 5 | `Homogeneous` | first-order **linear**; all additive terms have the same total degree in $(x,y)$ (degree-zero homogeneity), and a term free of $y$ is present |
| 6 | `LinearConstantCoeff` | order $\geq 2$; linear; the coefficients of the derivatives do not contain the independent variable |
| 7 | `CauchyEuler` | order $\geq 2$; the coefficient of $y^{(k)}$ is $c_k x^k$ |
| 8 | `ReductionOfOrder` | second-order linear, not of the types above |
| 9 | `PowerSeries` | fallback option for linear ODEs (first try a power series at an ordinary point, then Frobenius) |

`dsolve` tries the methods in this order; the first one that succeeds provides the result. The user can also specify a concrete method via the `hint` parameter.

### The Solver Architecture

```
dsolve(ctx, ode, hint?)
├── normalize_ode()          // simplify the equation
├── classify_ode()           // return the list of candidate methods
└── try in order:
    ├── first_order::solve_linear_first()
    ├── first_order::solve_bernoulli()
    ├── first_order::solve_separable()
    ├── first_order::solve_exact()
    ├── first_order::solve_homogeneous()
    ├── second_order::solve_constant_coeff()
    ├── second_order::solve_cauchy_euler()
    ├── second_order::solve_reduction_of_order()
    └── series::solve_power_series() or series::solve_frobenius()

dsolve_ivp(ctx, ode, y0, y1?)
└── laplace::solve_laplace()    // the Laplace transform method

dsolve_system(ctx, equations, funcs, var)
└── systems::solve_linear_system()  // 2×2 eigenvalue decomposition
```

### Return Value Types

The `ODESolution` enum covers all possible output forms:

| Variant | Meaning | Use case |
|---|---|---|
| `Explicit(expr)` | $y = \text{expr}$ | most successful solves |
| `Implicit(expr)` | $F(x,y) = 0$ | separable equations (implicit solutions) |
| `Parametric(x_t, y_t)` | $(x(t), y(t))$ | parametric solutions |
| `Series(expr, n_terms)` | truncated series | power series / Frobenius |
| `System([y1, y2])` | the individual components | $2 \times 2$ systems |
| `Unsolved(ode)` | the original equation | all methods failed |

### Integration with Symbolic Integration

The ODE solvers depend heavily on oCAS's symbolic integration capabilities:

- **First-order linear**: the integrating factor $\mu = \exp(\int p\, dx)$ and the final integral $\int \mu\, q\, dx$
- **Bernoulli**: the integrating factor of the linearized equation
- **Exact equations**: the double integral of the potential function
- **Homogeneous equations**: the integral after separation
- **Variation of parameters**: $\int y_i\, g / W\, dx$
- **Reduction of order**: $\int e^{-\int p\, dx} / y_1^2\, dx$

These integrals call oCAS's layered integration pipeline (rational function → Risch → trigonometric → special functions → heuristic); see [Symbolic Calculus](./symbolic-calculus.md) and [The Risch Integration Algorithm](./risch-algorithm.md).

### Limitations

Known limitations of the current implementation:

1. **Orders above two**: only linear constant-coefficient and Cauchy–Euler types are supported. General variable-coefficient higher-order equations have no dedicated solver.
2. **Nonlinear ODEs**: apart from Bernoulli (which linearizes), nonlinear equations usually return `Unsolved`.
3. **Laplace transform**: only handles constant-coefficient equations with integer coefficients, and the forcing term must have a tabulated transform.
4. **Series solutions**: expand 8 terms at $x_0 = 0$ by default. The Frobenius method only handles the case where the indicial equation has rational real roots.
5. **2×2 systems**: only constant coefficients are supported, and the eigenvalues must be integers or rationals.
6. **Verification**: `verify_solution` substitutes the solution back into the original equation and checks that the residual is zero (used only for testing).

---

## References

1. **Boyce, W. E. & DiPrima, R. C.** *Elementary Differential Equations and Boundary Value Problems*, 11th ed., Wiley, 2017. — Chapters 2–7 cover first-order equations (Ch.2), second-order linear (Ch.3–4), series solutions (Ch.5), Laplace transform (Ch.6), and systems (Ch.7).

2. **Coddington, E. A. & Levinson, N.** *Theory of Ordinary Differential Equations*, McGraw-Hill, 1955. — Rigorous treatment of existence, uniqueness, and stability.

3. **Ince, E. L.** *Ordinary Differential Equations*, Dover, 1956. — Classical reference on special functions and series methods.

4. **Zill, D. G.** *A First Course in Differential Equations with Modeling Applications*, 12th ed., Cengage, 2021. — Accessible introduction with applications.

5. **Tenenbaum, M. & Pollard, H.** *Ordinary Differential Equations*, Dover, 1985. — Comprehensive classical methods including exact equations and integrating factors.
