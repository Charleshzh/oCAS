# 数论 / Number Theory

`ocas_domain::number_theory` 模块提供计算数论栈：多模中国剩余定理、BPSW
素性判定、整数分解（试除、Pollard rho/p−1、Williams p+1、ECM）、离散对数
与经典积性函数。Python 与 C 绑定暴露同样的功能。

The `ocas_domain::number_theory` module provides the computational
number-theory stack: multi-modulus Chinese remaindering, BPSW primality,
integer factorization (trial division, Pollard rho/p−1, Williams p+1, ECM),
discrete logarithms, and the classical multiplicative functions. Python and C
bindings expose the same functionality.

## 素性判定 / Primality

- `is_prime(n)`：在 $n < 3.3\cdot10^{24}$ 内确定性成立。
- `primes::is_prime_bpsw(n)`：base-2 强 Miller–Rabin + 强 Lucas 测试
  （Selfridge 参数），目前无已知合数能通过。
- `primes::is_prime_u64(n)`：对整个 `u64` 范围确定性成立。

`is_prime(n)` is deterministic for $n < 3.3\cdot10^{24}$; `is_prime_bpsw`
combines a base-2 strong Miller–Rabin round with a strong Lucas test
(Selfridge parameters) — no composite is known to pass it; `is_prime_u64`
is deterministic over the whole `u64` range.

## 整数分解 / Integer Factorization

```rust
use ocas_domain::Integer;
use ocas_domain::number_theory::factor::factor_integer;

let f = factor_integer(&Integer::from(360));
// [(2, 3), (3, 2), (5, 1)]
```

`factor_integer` 先试除小因子，然后对每个合数余因子按递增光滑界升级
方法：先做一次快速 Pollard-rho 尝试，再逐轮执行 Pollard p−1 /
Williams p+1 / ECM（Suyama 参数化 + Montgomery 曲线 stage-1）。30 位
半素数在 release 模式下约一秒分解。

`factor_integer` peels small factors by trial division, then splits composite
cofactors with an escalating strategy: one quick Pollard-rho attempt, then
rounds of Pollard p−1 / Williams p+1 / ECM (Suyama parametrization,
Montgomery curves, stage 1) with growing smoothness bounds. A 30-digit
semiprime factors in about a second in release mode.

## 中国剩余定理 / Chinese Remainder Theorem

`crt::crt_many(&[(r1, m1), (r2, m2), ...])` 将任意同余式列表合并为单一
同余 `x ≡ R (mod M)`；模数无需两两互素，不一致的系统返回 `None`。

`crt::crt_many` merges any list of congruences into `x ≡ R (mod M)`; moduli
need not be coprime, and inconsistent systems return `None`.

## 离散对数 / Discrete Logarithms

- `dlog::dlog_bsgs(base, target, modulus)`：小步大步法，适用于小阶群。
- `dlog::dlog_pohlig_hellman(base, target, p)`：分解 `base` 在素数 `p`
  下的阶，对各素幂子群分别求解后经 CRT 合并。

`dlog_bsgs` is baby-step giant-step, practical for small groups;
`dlog_pohlig_hellman` factors the order of `base` modulo the prime `p` and
combines subgroup logarithms via CRT.

## 数论函数 / Number-Theoretic Functions

`functions::euler_phi`、`moebius_mu`、`divisor_tau`、`divisor_sigma(n, k)`、
`liouville_lambda`——均由 `|n|` 的素数分解式计算。二次剩余工具
（`legendre`、`jacobi`、Tonelli–Shanks 模平方根）位于父模块中。

All functions are computed from the prime factorization of `|n|`.
Quadratic-residue tools (`legendre`, `jacobi`, Tonelli–Shanks `mod_sqrt`)
live in the parent module.

## Python 绑定 / Python Bindings

```python
import ocas

ocas.factorint(360)            # [("2", 3), ("3", 2), ("5", 1)]
ocas.isprime(2**61 - 1)        # True
ocas.nextprime(10**6)          # 1000003
ocas.discrete_log(101, 2, 66)  # 83
ocas.crt([3, 5, 7], [2, 3, 2]) # (23, 105)
ocas.jacobi_symbol(2, 7)       # 1
ocas.totient(36)               # 12
ocas.mobius(30)                # -1
ocas.divisor_count(12)         # 6
ocas.divisor_sigma(12, 2)      # 210
ocas.liouville_lambda(12)      # -1
```

## C 绑定 / C Bindings

C API（`ocas_ntheory_*`，声明于 `include/ocas.h`）以十进制字符串传递
整数；返回的堆字符串用 `ocas_string_free` 释放。RAII 包装位于
`ocas::ntheory`（`include/ocas.hpp`）。

The C API (`ocas_ntheory_*`, declared in `include/ocas.h`) passes integers
as decimal strings; results are heap strings released with
`ocas_string_free`. RAII wrappers live in `ocas::ntheory`
(`include/ocas.hpp`).

## 模多项式 GCD / Modular Polynomial GCD

`ocas_poly::gcd::modular::gcd_modular_z` 用 Brown 模算法计算整系数稠密
单变量多项式的本原 GCD：对多个素数取 monic GCD 像，经 CRT（对称代表）
重构后用精确试除验证。它替代了朴素伪余式 GCD——后者在次数 ≳ 16 时
系数爆炸——并可处理 100 位系数的 50 次多项式。二元 `gcd_modular`
采用同样策略，另加主变量内容分离、monic 插值像与有理重构。

`gcd_modular_z` computes the primitive GCD of two dense univariate integer
polynomials with Brown's modular algorithm: monic GCD images modulo several
primes are combined with CRT (symmetric representatives) and confirmed by
exact trial division. It replaces the naive pseudo-remainder GCD, whose
coefficients explode for degrees ≳ 16, and handles degree-50 inputs with
100-digit coefficients. The bivariate `gcd_modular` applies the same strategy
with content separation, monic interpolation images, and rational
reconstruction.
