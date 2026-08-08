//! Rewrite rules for oCAS.
//!
//! A [`Rule`] pairs a pattern with a replacement builder. When the pattern
//! matches a sub-expression, the builder receives the wildcard bindings and
//! produces a replacement atom. Rules are applied by the [`simplify()`](crate::simplify::simplify)
//! function in a bottom-up traversal until no more changes occur.

use ocas_atom::{Atom, AtomArena, AtomNode, Symbol};
use ocas_core::FastHashMap;

use crate::matcher::{Bindings, MatchError, MatchValue, match_pattern};
use crate::pattern::{Pattern, PatternAlloc, WildcardLevel};

/// A rewrite rule.
///
/// Rules are typically created from parsed patterns via [`Rule::new`], or by
/// using the convenience constructors in the [`rules`](crate::rules) module.
///
/// The replacement and condition use boxed closures so that users can capture
/// their environment when building custom rules. Built-in rules still use
/// plain closures, so they incur only a single allocation per rule.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_atom::Symbol;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::matcher::{Bindings, MatchValue};
/// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
/// use ocas_rewrite::rules::Rule;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
/// let rule = Rule::new(pat, |bindings: &Bindings, _ctx: &AtomArena| {
///     let x = bindings.get(Symbol::new("x")).unwrap();
///     let MatchValue::Single(v) = x else { panic!("expected single"); };
///     _ctx.mul(&[_ctx.num(2), *v])
/// });
///
/// let y = ctx.var("y");
/// let result = rule.apply(&ctx, y).unwrap();
/// assert_eq!(result.to_string(), "2*y");
/// ```
pub struct Rule<'a> {
    pattern: Pattern<'a>,
    replacement: Replacement<'a>,
    condition: Option<Condition<'a>>,
    /// Optional template pattern; when present, `apply` instantiates it with
    /// the match bindings instead of calling `replacement`.
    template: Option<Pattern<'a>>,
}

type Replacement<'a> = Box<dyn Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a> + 'a>;
type Condition<'a> = Box<dyn Fn(&Bindings<'a>) -> bool + 'a>;

impl<'a> Rule<'a> {
    /// Create a rule from a pattern and a replacement builder.
    ///
    /// The `replacement` closure receives the bindings produced by a successful
    /// match and the arena context so it can construct new atoms.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_atom::Symbol;
    /// use ocas_core::arena::Arena;
    /// use ocas_rewrite::matcher::{Bindings, MatchValue};
    /// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
    /// use ocas_rewrite::rules::Rule;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
    /// let rule = Rule::new(pat, |bindings: &Bindings, _ctx: &AtomArena| {
    ///     let x = bindings.get(Symbol::new("x")).unwrap();
    ///     let MatchValue::Single(v) = x else { panic!("expected single"); };
    ///     _ctx.pow(*v, _ctx.num(2))
    /// });
    ///
    /// let z = ctx.var("z");
    /// let result = rule.apply(&ctx, z).unwrap();
    /// assert_eq!(result.to_string(), "z^2");
    /// ```
    pub fn new<F>(pattern: Pattern<'a>, replacement: F) -> Self
    where
        F: Fn(&Bindings<'a>, &AtomArena<'a>) -> Atom<'a> + 'a,
    {
        Self {
            pattern,
            replacement: Box::new(replacement),
            condition: None,
            template: None,
        }
    }

    /// Add a condition that must hold for the rule to fire.
    ///
    /// The condition is evaluated after a successful match but before the
    /// replacement is built.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_atom::AtomArena;
    /// use ocas_atom::Symbol;
    /// use ocas_core::arena::Arena;
    /// use ocas_rewrite::matcher::{Bindings, MatchValue};
    /// use ocas_rewrite::pattern::{Pattern, WildcardLevel};
    /// use ocas_rewrite::rules::Rule;
    ///
    /// let arena = Arena::new();
    /// let ctx = AtomArena::new(&arena);
    /// let pat = Pattern::Wildcard(Symbol::new("x"), WildcardLevel::Single);
    /// let rule = Rule::new(pat, |_bindings: &Bindings, ctx: &AtomArena| {
    ///     ctx.num(99)
    /// }).with_condition(|bindings: &Bindings| {
    ///     let x = bindings.get(Symbol::new("x")).unwrap();
    ///     let MatchValue::Single(v) = x else { return false; };
    ///     v.to_string() == "y"
    /// });
    ///
    /// let y = ctx.var("y");
    /// let z = ctx.var("z");
    /// assert_eq!(rule.apply(&ctx, y).unwrap().to_string(), "99");
    /// assert!(rule.apply(&ctx, z).is_none());
    /// ```
    pub fn with_condition<F>(mut self, condition: F) -> Self
    where
        F: Fn(&Bindings<'a>) -> bool + 'a,
    {
        self.condition = Some(Box::new(condition));
        self
    }

    /// Try to apply this rule to `atom`. Returns `Some` if the rule matched and
    /// the condition (if any) was satisfied.
    pub fn apply(&self, ctx: &'a AtomArena<'a>, atom: Atom<'a>) -> Option<Atom<'a>> {
        match match_pattern(self.pattern.clone(), atom) {
            Ok(bindings) => {
                if let Some(cond) = &self.condition
                    && !cond(&bindings)
                {
                    return None;
                }
                if let Some(tmpl) = &self.template {
                    let next = instantiate(ctx, tmpl, &bindings)?;
                    Some(ocas_atom::normalize::normalize(ctx, next))
                } else {
                    Some((self.replacement)(&bindings, ctx))
                }
            }
            Err(MatchError::NoMatch)
            | Err(MatchError::InconsistentBinding)
            | Err(MatchError::BudgetExhausted) => None,
        }
    }

    /// Build a rule from pattern/template strings (template instantiation).
    ///
    /// The template may contain the same wildcard names as the pattern;
    /// applying the rule instantiates the template with the match bindings and
    /// normalizes the result. This is the mechanism used by the integration
    /// rule library (ocas-calc) to express rules like
    /// `sin(x_)^n_ -> -sin(x_)^(n_-1)*cos(x_)/n_ + ...` without hand-written
    /// closures.
    pub fn from_template(
        ctx: &'a AtomArena<'a>,
        alloc: &'a impl PatternAlloc<'a>,
        pattern: &str,
        template: &str,
    ) -> Rule<'a> {
        let pat = pattern_from_str(ctx, alloc, pattern);
        let tmpl = pattern_from_str(ctx, alloc, template);
        Rule {
            pattern: pat,
            replacement: Box::new(|_, ctx| ctx.num(0)),
            condition: None,
            template: Some(tmpl),
        }
    }
}

/// Instantiate a template pattern with wildcard bindings.
///
/// Semantics:
/// - `Literal(a)` → `a`.
/// - `Wildcard(name, Single)` → the bound atom (unbound → `None`).
/// - `Wildcard(name, Sequence | NullSequence)` → only valid spliced into an
///   n-ary (`Add`/`Mul`/`Fun`) argument list, where the bound atoms are
///   inserted in order; at the top level → `None`.
/// - `Add`/`Mul`/`Fun` → rebuilt recursively (with sequence splice);
///   `Pow` → both sides rebuilt.
///
/// Returns `None` when a wildcard is unbound, bound at the wrong level, or
/// spliced outside an argument list.
pub fn instantiate<'a>(
    ctx: &'a AtomArena<'a>,
    template: &Pattern<'a>,
    bindings: &Bindings<'a>,
) -> Option<Atom<'a>> {
    match template {
        Pattern::Literal(a) => Some(*a),
        Pattern::Wildcard(name, WildcardLevel::Single) => match bindings.get(*name)? {
            MatchValue::Single(v) => Some(*v),
            MatchValue::Sequence(_) => None,
        },
        Pattern::Wildcard(name, level) => {
            // Sequence wildcards must be spliced inside an argument list;
            // a sequence wildcard at the root of the template is ambiguous.
            let _ = (name, level);
            None
        }
        Pattern::Add(pats) => {
            let args = instantiate_args(ctx, pats, bindings)?;
            Some(ctx.add(&args))
        }
        Pattern::Mul(pats) => {
            let args = instantiate_args(ctx, pats, bindings)?;
            Some(ctx.mul(&args))
        }
        Pattern::Fun(name, pats) => {
            let args = instantiate_args(ctx, pats, bindings)?;
            Some(ctx.fun(name.as_str(), &args))
        }
        Pattern::Pow(p_box) => {
            let (p_base, p_exp) = &**p_box;
            let base = instantiate(ctx, p_base, bindings)?;
            let exp = instantiate(ctx, p_exp, bindings)?;
            Some(ctx.pow(base, exp))
        }
    }
}

/// Instantiate a template argument list, splicing sequence wildcards.
fn instantiate_args<'a>(
    ctx: &'a AtomArena<'a>,
    pats: &[Pattern<'a>],
    bindings: &Bindings<'a>,
) -> Option<Vec<Atom<'a>>> {
    let mut out = Vec::with_capacity(pats.len());
    for pat in pats {
        match pat {
            Pattern::Wildcard(name, WildcardLevel::Sequence | WildcardLevel::NullSequence) => {
                match bindings.get(*name)? {
                    MatchValue::Sequence(slice) => out.extend_from_slice(slice),
                    MatchValue::Single(_) => return None,
                }
            }
            _ => out.push(instantiate(ctx, pat, bindings)?),
        }
    }
    Some(out)
}

/// Head key used to index rules for quick dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HeadKey {
    /// A function application with the given head symbol.
    Fun(Symbol),
    /// An addition (`Add`).
    Add,
    /// A product (`Mul`).
    Mul,
    /// A power (`Pow`).
    Pow,
    /// Anything else (numbers, variables, ...).
    Any,
}

/// The head key of an atom, used to look up candidate rules.
pub fn head_of(atom: Atom<'_>) -> HeadKey {
    match atom.node() {
        AtomNode::Fun(name, _) => HeadKey::Fun(*name),
        AtomNode::Add(_) => HeadKey::Add,
        AtomNode::Mul(_) => HeadKey::Mul,
        AtomNode::Pow(_, _) => HeadKey::Pow,
        AtomNode::Num(_) | AtomNode::Var(_) => HeadKey::Any,
    }
}

/// The bucket a rule's pattern root dispatches to.
fn pattern_head_key(pattern: &Pattern<'_>) -> HeadKey {
    match pattern {
        Pattern::Fun(name, _) => HeadKey::Fun(*name),
        Pattern::Add(_) => HeadKey::Add,
        Pattern::Mul(_) => HeadKey::Mul,
        Pattern::Pow(_) => HeadKey::Pow,
        Pattern::Literal(_) | Pattern::Wildcard(_, _) => HeadKey::Any,
    }
}

/// A head-indexed rule table.
///
/// [`RuleTable::from_rules`] buckets rules by their pattern's root node;
/// [`RuleTable::apply`] scans only the bucket matching the atom's head plus
/// the `Any` bucket. Within a bucket rules keep insertion order and the first
/// match wins — exactly the semantics of a linear scan over the same rules.
pub struct RuleTable<'a> {
    by_head: FastHashMap<HeadKey, Vec<Rule<'a>>>,
}

impl<'a> RuleTable<'a> {
    /// Build a table from rules, bucketing by pattern root.
    pub fn from_rules(rules: Vec<Rule<'a>>) -> Self {
        let mut by_head: FastHashMap<HeadKey, Vec<Rule<'a>>> = FastHashMap::default();
        for rule in rules {
            let key = pattern_head_key(&rule.pattern);
            by_head.entry(key).or_default().push(rule);
        }
        Self { by_head }
    }

    /// Apply the first matching rule to `atom` (head bucket, then `Any`).
    pub fn apply(&self, ctx: &'a AtomArena<'a>, atom: Atom<'a>) -> Option<Atom<'a>> {
        let key = head_of(atom);
        if let Some(rules) = self.by_head.get(&key) {
            for rule in rules {
                if let Some(next) = rule.apply(ctx, atom) {
                    return Some(next);
                }
            }
        }
        if key != HeadKey::Any
            && let Some(rules) = self.by_head.get(&HeadKey::Any)
        {
            for rule in rules {
                if let Some(next) = rule.apply(ctx, atom) {
                    return Some(next);
                }
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// Built-in algebraic rules
// ---------------------------------------------------------------------------

fn pattern_from_str<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
    s: &str,
) -> Pattern<'a> {
    let atom = ocas_parse::parse(ctx, s).expect("built-in rule pattern is valid");
    Pattern::from_atom(alloc, atom)
}

macro_rules! binding_single {
    ($bindings:expr, $name:expr) => {
        match $bindings.get(ocas_atom::Symbol::new($name)) {
            Some(crate::matcher::MatchValue::Single(atom)) => *atom,
            _ => panic!(concat!("expected single binding for '", $name, "'")),
        }
    };
}

/// `x + 0 -> x`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::add_zero;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.add(&[x, ctx.num(0)]);
/// let rule = add_zero(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "x");
/// ```
pub fn add_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ + 0"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `0 + x -> x`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::add_zero_left;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.add(&[ctx.num(0), x]);
/// let rule = add_zero_left(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "x");
/// ```
pub fn add_zero_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "0 + x_"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `x * 0 -> 0`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::mul_zero;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.mul(&[x, ctx.num(0)]);
/// let rule = mul_zero(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "0");
/// ```
pub fn mul_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ * 0"), |_bindings, ctx| {
        ctx.num(0)
    })
}

/// `0 * x -> 0`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::mul_zero_left;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.mul(&[ctx.num(0), x]);
/// let rule = mul_zero_left(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "0");
/// ```
pub fn mul_zero_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "0 * x_"), |_bindings, ctx| {
        ctx.num(0)
    })
}

/// `x * 1 -> x`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::mul_one;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.mul(&[x, ctx.num(1)]);
/// let rule = mul_one(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "x");
/// ```
pub fn mul_one<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ * 1"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `1 * x -> x`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::mul_one_left;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.mul(&[ctx.num(1), x]);
/// let rule = mul_one_left(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "x");
/// ```
pub fn mul_one_left<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "1 * x_"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// `x + x -> 2*x`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::add_same;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.add(&[x, x]);
/// let rule = add_same(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "2*x");
/// ```
pub fn add_same<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ + x_"), |bindings, ctx| {
        let x = binding_single!(bindings, "x");
        ctx.mul(&[ctx.num(2), x])
    })
}

/// `x ^ 0 -> 1`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::pow_zero;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.pow(x, ctx.num(0));
/// let rule = pow_zero(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "1");
/// ```
pub fn pow_zero<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ ^ 0"), |_bindings, ctx| {
        ctx.num(1)
    })
}

/// `x ^ 1 -> x`
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::pow_one;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let expr = ctx.pow(x, ctx.num(1));
/// let rule = pow_one(&ctx, &());
/// let result = rule.apply(&ctx, expr).unwrap();
/// assert_eq!(result.to_string(), "x");
/// ```
pub fn pow_one<'a>(ctx: &'a AtomArena<'a>, alloc: &'a impl PatternAlloc<'a>) -> Rule<'a> {
    Rule::new(pattern_from_str(ctx, alloc, "x_ ^ 1"), |bindings, _ctx| {
        binding_single!(bindings, "x")
    })
}

/// Return the default set of algebraic rewrite rules.
///
/// # Example
///
/// ```
/// use ocas_atom::AtomArena;
/// use ocas_core::arena::Arena;
/// use ocas_rewrite::rules::default_rules;
///
/// let arena = Arena::new();
/// let ctx = AtomArena::new(&arena);
/// let x = ctx.var("x");
/// let rules = default_rules(&ctx, &());
/// let expr = ctx.add(&[x, ctx.num(0)]);
/// let mut current = expr;
/// for rule in &rules {
///     if let Some(next) = rule.apply(&ctx, current) {
///         current = next;
///     }
/// }
/// assert_eq!(current.to_string(), "x");
/// ```
pub fn default_rules<'a>(
    ctx: &'a AtomArena<'a>,
    alloc: &'a impl PatternAlloc<'a>,
) -> Vec<Rule<'a>> {
    vec![
        add_zero(ctx, alloc),
        add_zero_left(ctx, alloc),
        mul_zero(ctx, alloc),
        mul_zero_left(ctx, alloc),
        mul_one(ctx, alloc),
        mul_one_left(ctx, alloc),
        add_same(ctx, alloc),
        pow_zero(ctx, alloc),
        pow_one(ctx, alloc),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    struct VecAlloc;

    impl<'a> PatternAlloc<'a> for VecAlloc {
        fn alloc_slice(&self, _items: &[Pattern<'a>]) -> &'a [Pattern<'a>] {
            // Not used by the matcher with the current Vec-based Pattern design.
            unreachable!()
        }
    }

    fn pat_atom<'a>(ctx: &'a AtomArena<'a>, alloc: &'a VecAlloc, s: &'a str) -> Pattern<'a> {
        let atom = ocas_parse::parse(ctx, s).unwrap();
        Pattern::from_atom(alloc, atom)
    }

    #[test]
    fn add_zero_applies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = add_zero(&ctx, &alloc);
        let x = ctx.var("x");
        let atom = ctx.add(&[x, ctx.num(0)]);
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result, x);
    }

    #[test]
    fn mul_zero_applies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = mul_zero(&ctx, &alloc);
        let x = ctx.var("x");
        let atom = ctx.mul(&[x, ctx.num(0)]);
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result, ctx.num(0));
    }

    #[test]
    fn add_same_applies() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = add_same(&ctx, &alloc);
        let x = ctx.var("x");
        let atom = ctx.add(&[x, x]);
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result.to_string(), "2*x");
    }

    #[test]
    fn rule_with_condition_can_block() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let pat = pat_atom(&ctx, &alloc, "x_");
        let rule = Rule::new(pat, |bindings, _ctx| binding_single!(bindings, "x"))
        .with_condition(|bindings| {
            matches!(
                bindings.get(ocas_atom::Symbol::new("x")),
                Some(crate::matcher::MatchValue::Single(a)) if matches!(a.node(), ocas_atom::AtomNode::Num(_))
            )
        });

        let x = ctx.var("x");
        assert!(rule.apply(&ctx, x).is_none());
        assert!(rule.apply(&ctx, ctx.num(5)).is_some());
    }

    #[test]
    fn from_template_instantiates_single_wildcards() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rule = Rule::from_template(&ctx, &alloc, "x_^n_", "x_^(n_+1)*(n_+1)^-1");
        let atom = ocas_parse::parse(&ctx, "x^3").unwrap();
        let result = rule.apply(&ctx, atom).unwrap();
        // n_ binds to 3, template is x^(n+1)*(n+1)^-1 = x^4/4.
        assert_eq!(result.to_string(), "(4^-1)*(x^4)");
        // Non-matching shape: no rewrite.
        assert!(rule.apply(&ctx, ctx.var("y")).is_none());
    }

    #[test]
    fn from_template_splices_sequence_wildcards() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        // `f(a__, b_)` -> `g(a__, b_)` with a non-empty sequence splice.
        let rule = Rule::from_template(&ctx, &alloc, "f(a__, b_)", "g(a__, b_)");
        let atom = ocas_parse::parse(&ctx, "f(x, y, z)").unwrap();
        let result = rule.apply(&ctx, atom).unwrap();
        assert_eq!(result.to_string(), "g(x, y, z)");
    }

    #[test]
    fn instantiate_unbound_wildcard_is_none() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let tmpl = pat_atom(&ctx, &alloc, "x_^2");
        let bindings = Bindings::new();
        assert!(instantiate(&ctx, &tmpl, &bindings).is_none());
    }

    #[test]
    fn rule_table_equals_linear_scan_on_default_rules() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let alloc = VecAlloc;
        let rules = default_rules(&ctx, &alloc);
        let table = RuleTable::from_rules(rules);
        let samples = [
            "x + 0", "0 + x", "x*0", "0*x", "x*1", "1*x", "x + x", "y^0", "y^1", "2*x + 0", "x*y",
            "0",
        ];
        for s in samples {
            let atom = ocas_parse::parse(&ctx, s).unwrap();
            let mut linear = atom;
            for rule in default_rules(&ctx, &alloc) {
                if let Some(next) = rule.apply(&ctx, linear) {
                    linear = next;
                }
            }
            let indexed = table.apply(&ctx, atom).unwrap_or(atom);
            assert_eq!(
                indexed, linear,
                "table vs linear disagree on {s}: {indexed} vs {linear}"
            );
        }
    }
}
