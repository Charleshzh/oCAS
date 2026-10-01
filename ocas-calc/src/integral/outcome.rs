//! Three-valued integration results and the certificate-gated entry point
//! (0.28.0).
//!
//! The classic [`integrate`](super::integrate) returns a single `Atom`: an
//! antiderivative when one was found, and the unevaluated `Integral(f, x)`
//! form otherwise. That loses information in two ways — it cannot say
//! whether a fallback is "I do not know" or "this is provably not
//! elementary", and it cannot say whether an emitted answer is *proved*
//! correct.
//!
//! [`integrate_outcome`] replaces the string criterion with a discipline:
//!
//! - [`Outcome::Found`] carries a machine-checkable
//!   [`Certificate`](super::certify::Certificate) — `D(F) − f` reduces to
//!   zero in the exact checker;
//! - [`Outcome::Unknown`] honestly returns the unevaluated form, and may
//!   additionally expose the pipeline's **uncertified** candidate;
//! - [`Outcome::ProvedNonElementary`] is reserved for the non-elementary
//!   layer (first producer planned with the `Li₂`/Meier-G wave); nothing in
//!   0.28.0 constructs it.
//!
//! The legacy `Atom`-returning entry points keep returning the pipeline's
//! candidate unchanged (including uncertified ones), so existing callers and
//! bindings keep their behaviour; they carry a rustdoc warning pointing at
//! [`integrate_outcome`] when a proof is required.

use ocas_atom::normalize::normalize;
use ocas_atom::{Atom, AtomArena, Symbol};

use super::certify::{CertDecline, Certificate, certify};
use super::{IntegrateOptions, integrate_with_options, residual};

/// Witness of a proof that an antiderivative does not exist in the
/// implemented elementary fragment.
///
/// Reserved: 0.28.0 defines the variant so the three-valued API is stable,
/// but no producer exists yet (the non-elementary wave adds the first one).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NonElementaryWitness<'a> {
    /// Name of the proof method that produced the witness.
    pub method: &'static str,
    /// The unevaluated `Integral(f, var)` form that is returned to callers.
    pub residue: Atom<'a>,
}

/// The result of a certified integration.
#[derive(Clone, Copy, Debug)]
pub enum Outcome<'a> {
    /// An antiderivative together with a machine-checkable certificate.
    Found {
        /// The antiderivative.
        value: Atom<'a>,
        /// The proof that `D(value) = f`.
        certificate: Certificate<'a>,
    },
    /// The integral was proved not to have an elementary antiderivative in
    /// the implemented fragment (no producer in 0.28.0).
    ProvedNonElementary {
        /// The proof witness, carrying the unevaluated form.
        witness: NonElementaryWitness<'a>,
    },
    /// Honest "I do not know": the unevaluated `Integral(f, var)` form.
    Unknown {
        /// The unevaluated residue, always an `Integral(f, var)` atom.
        residue: Atom<'a>,
        /// The pipeline's candidate, when it produced one but the exact
        /// checker could not certify it.
        ///
        /// **Not an answer.** It is surfaced for diagnostics, tracing and
        /// migration; using it without independent verification is exactly
        /// what this API exists to prevent.
        uncertified: Option<Atom<'a>>,
    },
}

impl<'a> Outcome<'a> {
    /// The value to show a caller: the certified antiderivative, or the
    /// honest residue.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::{AtomArena, Symbol};
    /// use ocas_calc::integral::outcome::integrate_outcome;
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let x = ctx.var("x");
    /// let two_x = ctx.mul(&[ctx.num(2), x]);
    /// let outcome = integrate_outcome(&ctx, two_x, Symbol::new("x"));
    /// assert!(outcome.is_found());
    /// assert!(!outcome.value().to_string().contains("Integral("));
    /// ```
    pub fn value(self) -> Atom<'a> {
        match self {
            Outcome::Found { value, .. } => value,
            Outcome::Unknown { residue, .. } => residue,
            Outcome::ProvedNonElementary { witness } => witness.residue,
        }
    }

    /// Whether the result carries a certificate.
    pub fn is_found(self) -> bool {
        matches!(self, Outcome::Found { .. })
    }

    /// The certificate, when there is one.
    pub fn certificate(self) -> Option<Certificate<'a>> {
        match self {
            Outcome::Found { certificate, .. } => Some(certificate),
            _ => None,
        }
    }

    /// The pipeline's uncertified candidate, when it produced one.
    pub fn uncertified(self) -> Option<Atom<'a>> {
        match self {
            Outcome::Unknown { uncertified, .. } => uncertified,
            _ => None,
        }
    }
}

/// Integrate `expr` and return a three-valued outcome, certifying the answer.
///
/// See the module docs for the emission discipline.
///
/// # Example
///
/// ```
/// use ocas_atom::{AtomArena, Symbol};
/// use ocas_calc::integral::outcome::{Outcome, integrate_outcome};
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// // ∫ exp(-x²) dx has no elementary antiderivative: the honest answer is
/// // the unevaluated form.
/// let e = ctx.fun("exp", &[ctx.mul(&[ctx.num(-1), ctx.pow(x, ctx.num(2))])]);
/// match integrate_outcome(&ctx, e, Symbol::new("x")) {
///     Outcome::Found { .. } => panic!("exp(-x²) must not be certified"),
///     Outcome::Unknown { residue, .. } => {
///         assert!(residue.to_string().contains("Integral("));
///     }
///     Outcome::ProvedNonElementary { .. } => panic!("no producer in 0.28.0"),
/// }
/// ```
pub fn integrate_outcome<'a>(ctx: &'a AtomArena<'a>, expr: Atom<'a>, var: Symbol) -> Outcome<'a> {
    integrate_outcome_with_options(ctx, expr, var, IntegrateOptions::default())
}

/// [`integrate_outcome`] with explicit pipeline options.
pub fn integrate_outcome_with_options<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    options: IntegrateOptions,
) -> Outcome<'a> {
    gated(ctx, expr, var, options).outcome
}

/// Diagnostic view of one certified integration: the raw pipeline candidate,
/// the exact checker's verdict, and the gated outcome.
#[doc(hidden)]
#[derive(Clone, Copy, Debug)]
pub struct DiagnosticOutcome<'a> {
    /// The pipeline's candidate, before any certification.
    pub raw: Atom<'a>,
    /// The exact checker's verdict for that candidate.
    pub certificate: Result<Certificate<'a>, CertDecline>,
    /// What [`integrate_outcome`] would return.
    pub outcome: Outcome<'a>,
}

/// Run the pipeline once and return both the raw candidate and the gated
/// outcome; used by the corpus harness to measure certificate coverage.
#[doc(hidden)]
pub fn integrate_diagnostic<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    options: IntegrateOptions,
) -> DiagnosticOutcome<'a> {
    gated(ctx, expr, var, options)
}

fn gated<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    var: Symbol,
    options: IntegrateOptions,
) -> DiagnosticOutcome<'a> {
    let normalized = normalize(ctx, expr);
    let raw = integrate_with_options(ctx, expr, var, options);
    let certificate = if residual::contains_residue(raw) {
        Err(CertDecline::NotInField)
    } else {
        certify(ctx, normalized, raw, var)
    };
    let outcome = match certificate {
        Ok(cert) => Outcome::Found {
            value: raw,
            certificate: cert,
        },
        Err(reason) => Outcome::Unknown {
            residue: if residual::contains_residue(raw) {
                raw
            } else {
                super::fallback(ctx, normalized, var)
            },
            uncertified: if residual::contains_residue(raw) {
                None
            } else {
                debug_assert!(matches!(
                    reason,
                    CertDecline::NotInField | CertDecline::Budget | CertDecline::NonZero
                ));
                Some(raw)
            },
        },
    };
    DiagnosticOutcome {
        raw,
        certificate,
        outcome,
    }
}

#[cfg(test)]
mod tests {
    use ocas_core::arena::Arena;
    use ocas_parse::parse;

    use super::*;

    #[test]
    fn certified_pairs_are_found() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        // ∫ 2x dx = x²
        let f = parse(&ctx, "2*x").expect("parse");
        let outcome = integrate_outcome(&ctx, f, Symbol::new("x"));
        assert!(outcome.is_found());
        assert!(outcome.certificate().is_some());
        assert!(!outcome.value().to_string().contains("Integral("));
    }

    #[test]
    fn non_elementary_integrals_are_unknown_with_a_residue() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let f = parse(&ctx, "exp(x^3)").expect("parse");
        let outcome = integrate_outcome(&ctx, f, Symbol::new("x"));
        assert!(!outcome.is_found());
        assert!(outcome.value().to_string().contains("Integral("));
    }

    #[test]
    fn diagnostic_reports_the_raw_candidate() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let f = parse(&ctx, "x^2*sin(x)").expect("parse");
        let d = integrate_diagnostic(&ctx, f, Symbol::new("x"), IntegrateOptions::default());
        // The raw candidate and the gated value agree for a certified case.
        assert!(!d.raw.to_string().contains("Integral("));
        assert_eq!(d.raw, d.outcome.value());
    }
}
