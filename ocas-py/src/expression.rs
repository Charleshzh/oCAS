//! Python `Expression` — a self-contained symbolic expression.
//!
//! Each [`Expression`] owns a private leaked `Arena` + `AtomArena<'static>`,
//! recovered on `Drop`. This mirrors the C API design and avoids cross-
//! reference lifetime entanglement between Python objects.

use ocas_atom::{Atom, AtomArena, Symbol, normalize::normalize};
use ocas_calc::{diff, integrate, substitute, taylor};
use ocas_core::arena::Arena;
use ocas_parse::parse;
use ocas_rewrite::rules::default_rules;
use ocas_rewrite::simplify::simplify;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Extend a string's lifetime to `'static`. Safe because oCAS atoms never
/// retain borrows of the input string — the parser copies characters into
/// arena-owned nodes.
///
/// # Safety
///
/// See above; only safe when the result is not stored beyond the input's
/// actual lifetime by code that depends on the borrow.
unsafe fn extend_str_lifetime(s: &str) -> &'static str {
    unsafe { std::mem::transmute::<&str, &'static str>(s) }
}

/// Internal storage behind an [`Expression`]: a leaked arena pair recovered
/// on drop.
struct ExprInner {
    arena_ptr: *mut Arena,
    ctx_ptr: *mut AtomArena<'static>,
    atom: Atom<'static>,
}

// SAFETY: the two heap allocations are not tied to any thread. The atom
// borrows them but they live until Drop.
unsafe impl Send for ExprInner {}
// SAFETY: pyo3 `#[pyclass]` (without `unordered`) requires `Send + Sync`.
// All `&self` method invocations are serialized by the GIL, so the
// `RefCell` inside `AtomArena` is never accessed concurrently.
// IMPORTANT: do not call `Python::allow_threads` with closures that access
// `ExprInner` — that would release the GIL and break this invariant.
unsafe impl Sync for ExprInner {}

impl Drop for ExprInner {
    fn drop(&mut self) {
        // SAFETY: both pointers came from `Box::into_raw`. Drop `ctx_ptr`
        // first because it borrows `arena_ptr`.
        unsafe {
            let _ = Box::from_raw(self.ctx_ptr);
            let _ = Box::from_raw(self.arena_ptr);
        }
    }
}

/// RAII guard that frees the leaked arena pair unless explicitly disarmed.
/// See [`ExprInner::build`] for usage.
struct ArenaGuard {
    arena_ptr: *mut Arena,
    ctx_ptr: *mut AtomArena<'static>,
    armed: bool,
}

impl ArenaGuard {
    fn new(arena_ptr: *mut Arena, ctx_ptr: *mut AtomArena<'static>) -> Self {
        ArenaGuard {
            arena_ptr,
            ctx_ptr,
            armed: true,
        }
    }

    fn disarm(&mut self) {
        self.armed = false;
    }
}

impl Drop for ArenaGuard {
    fn drop(&mut self) {
        if self.armed {
            // SAFETY: both pointers came from `Box::into_raw`.
            unsafe {
                let _ = Box::from_raw(self.ctx_ptr);
                let _ = Box::from_raw(self.arena_ptr);
            }
        }
    }
}

impl ExprInner {
    /// Borrow the atom arena as `&'static AtomArena<'static>`.
    fn ctx(&self) -> &'static AtomArena<'static> {
        // SAFETY: valid for as long as `ExprInner` is alive.
        unsafe { &*self.ctx_ptr }
    }

    /// Allocate a fresh arena pair.
    fn new_pair() -> (*mut Arena, *mut AtomArena<'static>) {
        let arena_box: Box<Arena> = Box::new(Arena::new());
        let arena_ptr = Box::into_raw(arena_box);
        // SAFETY: `arena_ptr` outlives `ExprInner`; recovered in Drop.
        let arena_ref: &'static Arena = unsafe { &*arena_ptr };
        let ctx = AtomArena::new(arena_ref);
        let ctx_ptr = Box::into_raw(Box::new(ctx));
        (arena_ptr, ctx_ptr)
    }

    /// Build from a closure that receives `&'static AtomArena<'static>`.
    fn build<F>(f: F) -> PyResult<Box<Self>>
    where
        F: FnOnce(&'static AtomArena<'static>) -> Result<Atom<'static>, String>,
    {
        let (arena_ptr, ctx_ptr) = Self::new_pair();
        let mut guard = ArenaGuard::new(arena_ptr, ctx_ptr);
        let ctx = unsafe { &*ctx_ptr };
        // If `f` or `normalize` panics, `guard` is dropped and frees the
        // arenas. On success we disarm and transfer ownership to ExprInner.
        let atom = f(ctx).map_err(PyValueError::new_err)?;
        let normalized = normalize(ctx, atom);
        guard.disarm();
        Ok(Box::new(ExprInner {
            arena_ptr,
            ctx_ptr,
            atom: normalized,
        }))
    }

    /// Parse a string.
    fn from_str(input: &str) -> PyResult<Box<Self>> {
        let static_input = unsafe { extend_str_lifetime(input) };
        Self::build(|ctx| match parse(ctx, static_input) {
            Ok(a) => Ok(a),
            Err(e) => Err(format!("parse error: {e}")),
        })
    }

    /// Rebuild from the string form of `src`.
    fn from_string_src(src: String) -> PyResult<Box<Self>> {
        let static_src = unsafe { extend_str_lifetime(&src) };
        Self::build(|ctx| match parse(ctx, static_src) {
            Ok(a) => Ok(a),
            Err(e) => Err(format!("parse error: {e}")),
        })
    }
}

/// A symbolic expression.
///
/// Construct from a string:
///
/// ```python
/// from ocas import Expression
/// e = Expression("x^2 + 2*x + 1")
/// print(e.diff("x"))
/// ```
#[pyclass(name = "Expression")]
pub struct Expression {
    inner: Box<ExprInner>,
}

#[pymethods]
impl Expression {
    /// Parse a string into an expression.
    #[new]
    fn new(input: &str) -> PyResult<Self> {
        Ok(Expression {
            inner: ExprInner::from_str(input)?,
        })
    }

    fn __str__(&self) -> String {
        self.inner.atom.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Expression({:?})", self.inner.atom.to_string())
    }

    fn __add__(&self, other: &Expression) -> PyResult<Expression> {
        let left = self.inner.atom.to_string();
        let right = other.inner.atom.to_string();
        let combined = format!("({left}) + ({right})");
        Ok(Expression {
            inner: ExprInner::from_string_src(combined)?,
        })
    }

    fn __sub__(&self, other: &Expression) -> PyResult<Expression> {
        let left = self.inner.atom.to_string();
        let right = other.inner.atom.to_string();
        let combined = format!("({left}) + (-1)*({right})");
        Ok(Expression {
            inner: ExprInner::from_string_src(combined)?,
        })
    }

    fn __mul__(&self, other: &Expression) -> PyResult<Expression> {
        let left = self.inner.atom.to_string();
        let right = other.inner.atom.to_string();
        let combined = format!("({left})*({right})");
        Ok(Expression {
            inner: ExprInner::from_string_src(combined)?,
        })
    }

    fn __pow__(&self, other: &Expression, _modulo: Option<&Expression>) -> PyResult<Expression> {
        let left = self.inner.atom.to_string();
        let right = other.inner.atom.to_string();
        let combined = format!("({left})^({right})");
        Ok(Expression {
            inner: ExprInner::from_string_src(combined)?,
        })
    }

    fn __neg__(&self) -> PyResult<Expression> {
        let src = self.inner.atom.to_string();
        Ok(Expression {
            inner: ExprInner::from_string_src(format!("(-1)*({src})"))?,
        })
    }

    fn __eq__(&self, other: &Expression) -> bool {
        // Compare normalized string forms.
        let a = normalize(self.inner.ctx(), self.inner.atom);
        let b = normalize(other.inner.ctx(), other.inner.atom);
        a.to_string() == b.to_string()
    }

    fn __hash__(&self) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut h = DefaultHasher::new();
        self.inner.atom.to_string().hash(&mut h);
        h.finish()
    }

    /// Return a copy of this expression.
    fn clone(&self) -> PyResult<Expression> {
        let src = self.inner.atom.to_string();
        Ok(Expression {
            inner: ExprInner::from_string_src(src)?,
        })
    }

    /// Simplify using the default rule set.
    fn simplify(&self) -> PyResult<Expression> {
        let src = self.inner.atom.to_string();
        let static_src = unsafe { extend_str_lifetime(&src) };
        ExprInner::build(|ctx| {
            let a = parse(ctx, static_src).map_err(|e| e.to_string())?;
            let rules = default_rules(ctx, &());
            Ok(simplify(ctx, a, &rules, 20))
        })
        .map(|inner| Expression { inner })
    }

    /// Differentiate with respect to `var`.
    fn diff(&self, var: &str) -> PyResult<Expression> {
        let src = self.inner.atom.to_string();
        let static_src = unsafe { extend_str_lifetime(&src) };
        let var_sym = Symbol::new(var);
        ExprInner::build(|ctx| match parse(ctx, static_src) {
            Ok(a) => Ok(diff(ctx, a, var_sym)),
            Err(e) => Err(e.to_string()),
        })
        .map(|inner| Expression { inner })
    }

    /// Integrate with respect to `var`.
    fn integrate(&self, var: &str) -> PyResult<Expression> {
        let src = self.inner.atom.to_string();
        let static_src = unsafe { extend_str_lifetime(&src) };
        let var_sym = Symbol::new(var);
        ExprInner::build(|ctx| match parse(ctx, static_src) {
            Ok(a) => Ok(integrate(ctx, a, var_sym)),
            Err(e) => Err(e.to_string()),
        })
        .map(|inner| Expression { inner })
    }

    /// Compute the Taylor series around `point` up to `order`.
    fn taylor(&self, var: &str, point: &Expression, order: usize) -> PyResult<Expression> {
        let expr_src = self.inner.atom.to_string();
        let point_src = point.inner.atom.to_string();
        let static_expr = unsafe { extend_str_lifetime(&expr_src) };
        let static_point = unsafe { extend_str_lifetime(&point_src) };
        let var_sym = Symbol::new(var);
        ExprInner::build(|ctx| {
            let e = parse(ctx, static_expr).map_err(|e| e.to_string())?;
            let p = parse(ctx, static_point).map_err(|e| e.to_string())?;
            Ok(taylor(ctx, e, var_sym, p, order))
        })
        .map(|inner| Expression { inner })
    }

    /// Substitute every occurrence of `var` with `replacement`.
    fn substitute(&self, var: &str, replacement: &Expression) -> PyResult<Expression> {
        let expr_src = self.inner.atom.to_string();
        let repl_src = replacement.inner.atom.to_string();
        let static_expr = unsafe { extend_str_lifetime(&expr_src) };
        let static_repl = unsafe { extend_str_lifetime(&repl_src) };
        let var_sym = Symbol::new(var);
        ExprInner::build(|ctx| {
            let e = parse(ctx, static_expr).map_err(|e| e.to_string())?;
            let r = parse(ctx, static_repl).map_err(|e| e.to_string())?;
            Ok(substitute(ctx, e, var_sym, r))
        })
        .map(|inner| Expression { inner })
    }
}
