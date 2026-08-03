# Number Theory

The `ocas_domain::number_theory` module provides the computational
number-theory stack: multi-modulus Chinese remaindering, BPSW primality,
integer factorization (trial division, Pollard rho/p−1, Williams p+1, ECM),
discrete logarithms, and the classical multiplicative functions. Python and C
bindings expose the same functionality.

## Primality Testing

- `is_prime(n)`: deterministic for $n < 3.3\cdot10^{24}$.
- `primes::is_prime_bpsw(n)`: base-2 strong Miller–Rabin + strong Lucas test
  (Selfridge parameters); no composite is known to pass it.
- `primes::is_prime_u64(n)`: deterministic over the whole `u64` range.

## Integer Factorization

```rust
use ocas_domain::Integer;
use ocas_domain::number_theory::factor::factor_integer;

let f = factor_integer(&Integer::from(360));
// [(2, 3), (3, 2), (5, 1)]
```

`factor_integer` peels small factors by trial division (primes up to 1000),
then splits each composite cofactor with an escalating strategy: one quick
Pollard-rho attempt, then rounds of Pollard p−1 / Williams p+1 / ECM
(Suyama parametrization, Montgomery curves, stage 1) with growing smoothness
bounds (the curve budget scales as ≈ B1/550). A 30-digit semiprime factors
in about a second in release mode.

## Chinese Remainder Theorem

`crt::crt_many(&[(r1, m1), (r2, m2), ...])` merges any list of congruences
into a single congruence `x ≡ R (mod M)`; the moduli need not be pairwise
coprime, and inconsistent systems return `None`.

## Discrete Logarithms

- `dlog::dlog_bsgs(base, target, modulus)`: baby-step giant-step, practical
  for small groups.
- `dlog::dlog_pohlig_hellman(base, target, p)`: factors the order of `base`
  modulo the prime `p` and combines the subgroup logarithms via CRT.

## Number-Theoretic Functions

`functions::euler_phi`, `moebius_mu`, `divisor_tau`, `divisor_sigma(n, k)`,
`liouville_lambda` — all computed from the prime factorization of `|n|`.
Quadratic-residue tools (`legendre`, `jacobi`, Tonelli–Shanks `mod_sqrt`)
live in the parent module.

## Python Bindings

Note the argument order of two functions: `discrete_log(p, base, target)`
takes the modulus first (it solves `base^x ≡ target (mod p)` for `x`);
`crt(moduli, residues)` takes the modulus list first and the residue list
second.

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

## C Bindings

The C API (`ocas_ntheory_*`, declared in `include/ocas.h`) passes integers as
decimal strings; results are heap strings released with `ocas_string_free`.
RAII wrappers live in `ocas::ntheory` (`include/ocas.hpp`).

## Modular Polynomial GCD

`ocas_poly::gcd::modular::gcd_modular_z` computes the primitive GCD of two
dense univariate integer polynomials with Brown's modular algorithm: monic GCD
images modulo several primes are combined with CRT (symmetric representatives)
and confirmed by exact trial division. It replaces the naive pseudo-remainder
GCD, whose coefficients explode for degrees ≳ 16, and handles degree-50 inputs
with 100-digit coefficients. The bivariate `gcd_modular` applies the same
strategy with content separation, monic interpolation images, and rational
reconstruction.
