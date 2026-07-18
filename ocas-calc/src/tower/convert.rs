//! Bidirectional conversion between [`Atom`] expressions and multivariate
//! rational functions over `ℚ`.
//!
//! The Risch integrator works on elements of a differential field tower
//! `ℚ(x, t₁, …, tₙ)`. Elements are represented as
//! [`RationalPolynomial<Rational, Lex>`] whose variables are the tower
//! generators. This module converts between the atom world (what the
//! parser and printer understand) and the polynomial world (what the
//! algorithms operate on).
//!
//! Generators are arbitrary atoms — typically the integration variable
//! `x` and function applications such as `log(x)` or `exp(x)`. They are
//! listed from the bottom to the top of the tower; the *last* generator
//! is the main variable for univariate views.
//!
//! Rational coefficients follow the atom convention `p * q^-1`: a
//! denominator is represented as a negative power, matching the rest of
//! the calculus modules.

use ocas_atom::{Atom, AtomArena, AtomNode};
use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::{Lex, RationalPolynomial, SparseMultivariatePolynomial};

/// A rational function over `ℚ` in the tower generators.
pub type GeneratorField = RationalPolynomial<RationalDomain, Lex>;

type Poly = SparseMultivariatePolynomial<RationalDomain, Lex>;

/// Convert `atom` into a rational function over the given generators.
///
/// Returns `None` when the atom is not a rational function of the
/// generators: an unknown variable or function application, a
/// non-integer exponent, or a negative power of zero.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_calc::tower::convert::atom_to_rational;
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.num(1)]);
/// let rf = atom_to_rational(expr, &[x]).unwrap();
/// assert_eq!(rf.numerator.n_terms(), 2);
/// ```
pub fn atom_to_rational<'a>(atom: Atom<'a>, gens: &[Atom<'a>]) -> Option<GeneratorField> {
    atom_to_rational_extended(atom, gens, gens.len())
}

/// Like [`atom_to_rational`] but embeds into a polynomial ring with
/// `n_vars >= gens.len()` variables (extra variables are left unused).
///
/// Used by the tower construction, where generator arguments are
/// converted before the full set of generators is known.
pub(crate) fn atom_to_rational_extended<'a>(
    atom: Atom<'a>,
    gens: &[Atom<'a>],
    n_vars: usize,
) -> Option<GeneratorField> {
    debug_assert!(n_vars >= gens.len());
    match atom.node() {
        AtomNode::Num(n) => Some(constant(*n, n_vars)),
        AtomNode::Var(_) | AtomNode::Fun(..) => {
            let idx = gens.iter().position(|g| *g == atom)?;
            Some(variable(idx, n_vars))
        }
        AtomNode::Add(args) => {
            let mut acc = GeneratorField::zero(&RationalDomain, n_vars);
            for a in args.iter() {
                acc = acc.add(&atom_to_rational_extended(*a, gens, n_vars)?);
            }
            Some(acc)
        }
        AtomNode::Mul(args) => {
            let mut acc = GeneratorField::one(&RationalDomain, n_vars);
            for a in args.iter() {
                acc = acc.mul(&atom_to_rational_extended(*a, gens, n_vars)?);
            }
            Some(acc)
        }
        AtomNode::Pow(base, exp) => {
            let AtomNode::Num(n) = exp.node() else {
                return None;
            };
            let b = atom_to_rational_extended(*base, gens, n_vars)?;
            if *n >= 0 {
                Some(pow_u64(&b, *n as u64))
            } else {
                GeneratorField::one(&RationalDomain, n_vars).div(&pow_u64(&b, n.unsigned_abs()))
            }
        }
    }
}

/// Convert a rational function back into an atom over the generators.
///
/// Terms are emitted in deterministic order: the last generator is the
/// most significant, and exponents descend. Returns `None` when a
/// coefficient does not fit into an `i64`.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_calc::tower::convert::{atom_to_rational, rational_to_atom};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.add(&[ctx.pow(x, ctx.num(2)), ctx.mul(&[ctx.num(2), x]), ctx.num(1)]);
/// let rf = atom_to_rational(expr, &[x]).unwrap();
/// let back = rational_to_atom(&ctx, &rf, &[x]).unwrap();
/// assert_eq!(back.to_string(), "(x^2) + (2*x) + 1");
/// ```
pub fn rational_to_atom<'a>(
    ctx: &'a AtomArena<'a>,
    rf: &GeneratorField,
    gens: &[Atom<'a>],
) -> Option<Atom<'a>> {
    let num = poly_to_atom(ctx, &rf.numerator, gens)?;
    if is_const_one(&rf.denominator) {
        return Some(num);
    }
    // 1/den prints as den^-1 without a leading 1* factor.
    if is_const_one(&rf.numerator) {
        let den = poly_to_atom(ctx, &rf.denominator, gens)?;
        return Some(ctx.pow(den, ctx.num(-1)));
    }
    let den = poly_to_atom(ctx, &rf.denominator, gens)?;
    Some(ctx.mul(&[num, ctx.pow(den, ctx.num(-1))]))
}

fn constant(n: i64, n_vars: usize) -> GeneratorField {
    if n == 0 {
        return GeneratorField::zero(&RationalDomain, n_vars);
    }
    let poly = Poly::from_terms(
        RationalDomain,
        n_vars,
        vec![(vec![0; n_vars], Rational::new(n, 1))],
    );
    GeneratorField::from_polynomial(poly)
}

fn variable(idx: usize, n_vars: usize) -> GeneratorField {
    let mut exp = vec![0usize; n_vars];
    exp[idx] = 1;
    let poly = Poly::from_terms(RationalDomain, n_vars, vec![(exp, RationalDomain.one())]);
    GeneratorField::from_polynomial(poly)
}

fn pow_u64(base: &GeneratorField, mut exp: u64) -> GeneratorField {
    let mut result = GeneratorField::one(&RationalDomain, base.n_vars());
    let mut b = base.clone();
    while exp > 0 {
        if exp & 1 == 1 {
            result = result.mul(&b);
        }
        b = b.mul(&b);
        exp >>= 1;
    }
    result
}

fn is_const_one(p: &Poly) -> bool {
    p.n_terms() == 1 && p.domain().is_one(&p.coeff(&vec![0; p.n_vars()]))
}

/// Convert a rational constant to an atom using the `p * q^-1` convention.
///
/// Returns `None` when the numerator or denominator does not fit into an
/// `i64`.
pub fn rational_const_to_atom<'a>(ctx: &'a AtomArena<'a>, coeff: &Rational) -> Option<Atom<'a>> {
    term_to_atom(ctx, &[], coeff, &[])
}

fn poly_to_atom<'a>(ctx: &'a AtomArena<'a>, poly: &Poly, gens: &[Atom<'a>]) -> Option<Atom<'a>> {
    let mut terms: Vec<_> = poly.terms_ref().iter().collect();
    // Deterministic order: the last generator is most significant,
    // exponents in descending order.
    terms.sort_by(|(e1, _), (e2, _)| e2.iter().rev().cmp(e1.iter().rev()));
    let mut sum = Vec::with_capacity(terms.len());
    for (exp, coeff) in terms {
        sum.push(term_to_atom(ctx, exp, coeff, gens)?);
    }
    Some(match sum.len() {
        0 => ctx.num(0),
        1 => sum[0],
        _ => ctx.add(&sum),
    })
}

fn term_to_atom<'a>(
    ctx: &'a AtomArena<'a>,
    exp: &[usize],
    coeff: &Rational,
    gens: &[Atom<'a>],
) -> Option<Atom<'a>> {
    let p = coeff.numer().to_i64()?;
    let q = coeff.denom().to_i64()?;
    let mut factors: Vec<Atom> = Vec::new();
    let has_monomial = exp.iter().any(|&e| e > 0);
    if q != 1 {
        if p != 1 {
            factors.push(ctx.num(p));
        }
        factors.push(ctx.pow(ctx.num(q), ctx.num(-1)));
    } else if p != 1 || !has_monomial {
        factors.push(ctx.num(p));
    }
    for (i, &e) in exp.iter().enumerate() {
        if e == 0 {
            continue;
        }
        let g = gens[i];
        factors.push(if e == 1 {
            g
        } else {
            ctx.pow(g, ctx.num(e as i64))
        });
    }
    Some(match factors.len() {
        0 => ctx.num(1),
        1 => factors[0],
        _ => ctx.mul(&factors),
    })
}

#[cfg(test)]
mod tests {
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    use super::*;

    #[test]
    fn polynomial_roundtrip() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let expr = ctx.add(&[
            ctx.pow(x, ctx.num(2)),
            ctx.mul(&[ctx.num(2), x]),
            ctx.num(1),
        ]);
        let rf = atom_to_rational(expr, &[x]).unwrap();
        assert!(is_const_one(&rf.denominator));
        let back = rational_to_atom(&ctx, &rf, &[x]).unwrap();
        assert_eq!(back.to_string(), "(x^2) + (2*x) + 1");
    }

    #[test]
    fn rational_coefficient_roundtrip() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // (1/2)*x + 3*x^-1
        let half_x = ctx.mul(&[ctx.pow(ctx.num(2), ctx.num(-1)), x]);
        let three_over_x = ctx.mul(&[ctx.num(3), ctx.pow(x, ctx.num(-1))]);
        let rf = atom_to_rational(ctx.add(&[half_x, three_over_x]), &[x]).unwrap();
        // Canonical form: ((1/2)*x^2 + 3) / x (monic denominator).
        let back = rational_to_atom(&ctx, &rf, &[x]).unwrap();
        assert_eq!(back.to_string(), "(((2^-1)*(x^2)) + 3)*(x^-1)");
        // Atom -> rf -> atom -> rf is a fixed point.
        let rf2 = atom_to_rational(back, &[x]).unwrap();
        assert_eq!(rf, rf2);
    }

    #[test]
    fn function_generators() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let log_x = ctx.fun("log", &[x]);
        // x*log(x) + 1 over generators [x, log(x)]
        let expr = ctx.add(&[ctx.mul(&[x, log_x]), ctx.num(1)]);
        let rf = atom_to_rational(expr, &[x, log_x]).unwrap();
        let back = rational_to_atom(&ctx, &rf, &[x, log_x]).unwrap();
        assert_eq!(back.to_string(), "(x*(log(x))) + 1");
    }

    #[test]
    fn negative_exponent_denominator() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // (x+1)^-2 -> numerator 1, denominator (x+1)^2
        let expr = ctx.pow(ctx.add(&[x, ctx.num(1)]), ctx.num(-2));
        let rf = atom_to_rational(expr, &[x]).unwrap();
        assert!(is_const_one(&rf.numerator));
        let back = rational_to_atom(&ctx, &rf, &[x]).unwrap();
        assert_eq!(back.to_string(), "((x^2) + (2*x) + 1)^-1");
    }

    #[test]
    fn rejects_non_rational_input() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        // Unknown variable.
        assert!(atom_to_rational(ctx.var("y"), &[x]).is_none());
        // Unknown function application.
        assert!(atom_to_rational(ctx.fun("sin", &[x]), &[x]).is_none());
        // Non-integer exponent.
        assert!(atom_to_rational(ctx.pow(x, ctx.var("y")), &[x]).is_none());
        // Negative power of zero is undefined.
        assert!(atom_to_rational(ctx.pow(ctx.num(0), ctx.num(-1)), &[x]).is_none());
    }

    #[test]
    fn zero_and_constants() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let x = ctx.var("x");
        let zero = atom_to_rational(ctx.num(0), &[x]).unwrap();
        assert!(zero.is_zero());
        assert_eq!(
            rational_to_atom(&ctx, &zero, &[x]).unwrap().to_string(),
            "0"
        );
        let five = atom_to_rational(ctx.num(5), &[x]).unwrap();
        assert_eq!(
            rational_to_atom(&ctx, &five, &[x]).unwrap().to_string(),
            "5"
        );
    }
}
