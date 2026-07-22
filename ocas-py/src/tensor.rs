//! Python bindings for the basic tensor algebra module.
//!
//! Wraps [`ocas_atom::tensor`] — an independent `Tensor` type with index
//! slots, variance, and slot symmetry, plus explicit contraction and a
//! symmetrisation sign. Each [`Tensor`][PyTensor] owns a private leaked
//! arena pair (mirroring [`Expression`](crate::expression::Expression)).
//!
//! ```python
//! from ocas import Tensor, contract_tensors, tensor_symmetrise_sign
//!
//! # T^i_j · U^j_k = (TU)^i_k  (partial contraction over j)
//! t = Tensor("T", [("i", "upper"), ("j", "lower")])
//! u = Tensor("U", [("j", "upper"), ("k", "lower")])
//! kind, payload = contract_tensors(t, u)
//! assert kind == "product"
//!
//! # Antisymmetric ε_ab has a sign under slot swap.
//! eps = Tensor("eps", [("a", "lower"), ("b", "lower")], symmetry="antisymmetric")
//! assert tensor_symmetrise_sign(eps) in (1, -1)
//! ```

use ocas_atom::tensor::{
    Contracted, IndexPosition, IndexSlot, Symmetry, Tensor, contract, symmetrise_sign,
};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

// ------------------------------------------------------------------
//  Arena management (leaked pair, recovered on Drop)
// ------------------------------------------------------------------

/// Extend a string's lifetime to `'static`. Safe because atoms never retain
/// borrows of the input string.
unsafe fn extend_str_lifetime(s: &str) -> &'static str {
    unsafe { std::mem::transmute::<&str, &'static str>(s) }
}

/// Internal storage behind a [`PyTensor`]: a leaked arena pair recovered on
/// drop. See [`crate::expression::ExprInner`] for the same pattern.
struct TensorInner {
    arena_ptr: *mut Arena,
    ctx_ptr: *mut AtomArena<'static>,
    tensor: Tensor<'static>,
}

// SAFETY: matches crate::expression::ExprInner — the heap allocations are not
// tied to any thread, and pyclass method invocations are GIL-serialized.
unsafe impl Send for TensorInner {}
unsafe impl Sync for TensorInner {}

impl Drop for TensorInner {
    fn drop(&mut self) {
        // SAFETY: both pointers came from `Box::into_raw`. Drop `ctx_ptr`
        // first because it borrows `arena_ptr`.
        unsafe {
            let _ = Box::from_raw(self.ctx_ptr);
            let _ = Box::from_raw(self.arena_ptr);
        }
    }
}

impl TensorInner {
    /// Borrow the atom arena as `&'static AtomArena<'static>`.
    fn ctx(&self) -> &'static AtomArena<'static> {
        // SAFETY: valid while `TensorInner` is alive.
        unsafe { &*self.ctx_ptr }
    }

    /// Build a fresh arena pair.
    fn new_pair() -> (*mut Arena, *mut AtomArena<'static>) {
        let arena_box: Box<Arena> = Box::new(Arena::new());
        let arena_ptr = Box::into_raw(arena_box);
        // SAFETY: `arena_ptr` outlives `TensorInner`; recovered in Drop.
        let arena_ref: &'static Arena = unsafe { &*arena_ptr };
        let ctx = AtomArena::new(arena_ref);
        let ctx_ptr = Box::into_raw(Box::new(ctx));
        (arena_ptr, ctx_ptr)
    }

    /// Build a `TensorInner` from a closure that constructs the tensor in the
    /// freshly-allocated arena.
    fn build<F>(f: F) -> PyResult<Box<Self>>
    where
        F: FnOnce(&'static AtomArena<'static>) -> Tensor<'static>,
    {
        let (arena_ptr, ctx_ptr) = Self::new_pair();
        // If `f` panics we must free the arenas. Use a guard.
        struct Guard {
            arena_ptr: *mut Arena,
            ctx_ptr: *mut AtomArena<'static>,
            armed: bool,
        }
        impl Drop for Guard {
            fn drop(&mut self) {
                if self.armed {
                    unsafe {
                        let _ = Box::from_raw(self.ctx_ptr);
                        let _ = Box::from_raw(self.arena_ptr);
                    }
                }
            }
        }
        let mut g = Guard {
            arena_ptr,
            ctx_ptr,
            armed: true,
        };
        let ctx = unsafe { &*ctx_ptr };
        let tensor = f(ctx);
        g.armed = false;
        Ok(Box::new(TensorInner {
            arena_ptr,
            ctx_ptr,
            tensor,
        }))
    }
}

// ------------------------------------------------------------------
//  Helpers
// ------------------------------------------------------------------

/// Parse `"upper"` / `"lower"` into an [`IndexPosition`].
fn parse_position(s: &str) -> PyResult<IndexPosition> {
    match s.to_ascii_lowercase().as_str() {
        "upper" | "up" | "contravariant" => Ok(IndexPosition::Upper),
        "lower" | "down" | "covariant" => Ok(IndexPosition::Lower),
        _ => Err(PyValueError::new_err(format!(
            "position must be 'upper' or 'lower', got {s:?}"
        ))),
    }
}

fn position_str(p: IndexPosition) -> &'static str {
    match p {
        IndexPosition::Upper => "upper",
        IndexPosition::Lower => "lower",
    }
}

fn parse_symmetry(s: &str) -> PyResult<Symmetry> {
    match s.to_ascii_lowercase().as_str() {
        "none" | "" => Ok(Symmetry::None),
        "symmetric" | "sym" => Ok(Symmetry::Symmetric),
        "antisymmetric" | "antisym" | "skew" => Ok(Symmetry::Antisymmetric),
        _ => Err(PyValueError::new_err(format!(
            "symmetry must be 'none', 'symmetric', or 'antisymmetric', got {s:?}"
        ))),
    }
}

fn symmetry_str(s: Symmetry) -> &'static str {
    match s {
        Symmetry::None => "none",
        Symmetry::Symmetric => "symmetric",
        Symmetry::Antisymmetric => "antisymmetric",
    }
}

/// A tensor: a named object with a list of index slots and an optional slot
/// symmetry. Each tensor owns its own private arena; contraction rebuilds the
/// operands into a fresh arena so lifetimes stay decoupled.
#[pyclass(name = "Tensor")]
pub struct PyTensor {
    inner: Box<TensorInner>,
}

#[pymethods]
impl PyTensor {
    /// Create a tensor from a name and a list of `(label, position)` slots.
    ///
    /// `position` is the string `"upper"` (or `"up"`, `"contravariant"`) or
    /// `"lower"` (or `"down"`, `"covariant"`). The optional `symmetry`
    /// keyword is one of `"none"` (default), `"symmetric"`, or
    /// `"antisymmetric"`.
    #[new]
    #[pyo3(signature = (name, slots, symmetry="none"))]
    fn new(name: &str, slots: &Bound<'_, PyAny>, symmetry: &str) -> PyResult<Self> {
        let sym = parse_symmetry(symmetry)?;
        let parsed: Vec<(String, IndexPosition)> = slots
            .try_iter()
            .map_err(|_| PyValueError::new_err("slots must be a list of (label, position) pairs"))?
            .map(|item| -> PyResult<(String, IndexPosition)> {
                let item = item?;
                let (label, pos): (String, String) = item.extract().map_err(|_| {
                    PyValueError::new_err("each slot must be a (label, position) pair")
                })?;
                Ok((label, parse_position(&pos)?))
            })
            .collect::<PyResult<_>>()?;
        let static_name = unsafe { extend_str_lifetime(name) };
        let symbol = Symbol::new(static_name);
        let inner = TensorInner::build(|ctx| {
            let slots: Vec<IndexSlot<'static>> = parsed
                .iter()
                .map(|(label, pos)| {
                    let static_label = unsafe { extend_str_lifetime(label) };
                    IndexSlot::new(ctx.var(static_label), *pos)
                })
                .collect();
            Tensor::new(symbol, slots).with_symmetry(sym)
        })?;
        Ok(PyTensor { inner })
    }

    /// The tensor name.
    #[getter]
    fn name(&self) -> String {
        self.inner.tensor.name().as_str().to_string()
    }

    /// The tensor arity (number of slots).
    #[getter]
    fn rank(&self) -> usize {
        self.inner.tensor.rank()
    }

    /// The slot symmetry string ("none", "symmetric", or "antisymmetric").
    #[getter]
    fn symmetry(&self) -> &'static str {
        symmetry_str(self.inner.tensor.symmetry())
    }

    /// Return the slots as a list of `(label, position)` string pairs.
    fn slots(&self) -> Vec<(String, &'static str)> {
        self.inner
            .tensor
            .slots()
            .iter()
            .map(|s| (s.label().to_string(), position_str(s.position())))
            .collect()
    }

    /// Return the dummy labels (labels occurring exactly twice across the
    /// slots).
    fn dummy_labels(&self) -> Vec<String> {
        self.inner
            .tensor
            .dummy_labels()
            .into_iter()
            .map(|a| a.to_string())
            .collect()
    }

    /// Render the tensor as an `Atom` function node `name(slot, slot, ...)`.
    fn to_string_atom(&self) -> String {
        self.inner.tensor.to_atom(self.inner.ctx()).to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "Tensor({:?}, rank={}, symmetry={:?})",
            self.name(),
            self.rank(),
            self.symmetry()
        )
    }
}

// ------------------------------------------------------------------
//  contract and symmetrise_sign
// ------------------------------------------------------------------

/// Build an independent [`PyTensor`] (own arena) from name/slots/symmetry.
fn rebuild_tensor(
    name: &str,
    sym: Symmetry,
    slots: &[(String, IndexPosition)],
) -> PyResult<PyTensor> {
    let static_name = unsafe { extend_str_lifetime(name) };
    let inner = TensorInner::build(|ctx| {
        let slots: Vec<IndexSlot<'static>> = slots
            .iter()
            .map(|(label, pos)| {
                let static_label = unsafe { extend_str_lifetime(label) };
                IndexSlot::new(ctx.var(static_label), *pos)
            })
            .collect();
        Tensor::new(Symbol::new(static_name), slots).with_symmetry(sym)
    })?;
    Ok(PyTensor { inner })
}

/// Snapshot a tensor's name, symmetry, and slots as plain `String`/enum data
/// so it can be rebuilt into a fresh arena.
fn snapshot(tensor: &Tensor<'_>) -> (String, Symmetry, Vec<(String, IndexPosition)>) {
    let name = tensor.name().as_str().to_string();
    let sym = tensor.symmetry();
    let slots: Vec<(String, IndexPosition)> = tensor
        .slots()
        .iter()
        .map(|s| (s.label().to_string(), s.position()))
        .collect();
    (name, sym, slots)
}

/// Contract two tensors by summing over shared dummy indices (equal label,
/// opposite variance).
///
/// Returns a `(kind, payload)` tuple where `kind` is `"product"` or
/// `"scalar"`. For `"product"`, `payload` is a list of resulting tensors
/// (their free slots concatenated). For `"scalar"`, `payload` is the
/// string form of the contracted atom expression.
#[pyfunction]
pub fn contract_tensors<'py>(
    py: Python<'py>,
    a: &PyTensor,
    b: &PyTensor,
) -> PyResult<Bound<'py, PyAny>> {
    // Allocate a single shared arena for the contraction computation. It is
    // dropped before this function returns; the result PyTensors are rebuilt
    // into independent arenas via `rebuild_tensor`.
    let (arena_ptr, ctx_ptr) = TensorInner::new_pair();
    struct DropGuard {
        arena_ptr: *mut Arena,
        ctx_ptr: *mut AtomArena<'static>,
    }
    impl Drop for DropGuard {
        fn drop(&mut self) {
            unsafe {
                let _ = Box::from_raw(self.ctx_ptr);
                let _ = Box::from_raw(self.arena_ptr);
            }
        }
    }
    let _guard = DropGuard { arena_ptr, ctx_ptr };
    let ctx: &'static AtomArena<'static> = unsafe { &*ctx_ptr };

    // Rebuild a and b into the shared arena.
    let (a_name, a_sym, a_slots_data) = snapshot(&a.inner.tensor);
    let (b_name, b_sym, b_slots_data) = snapshot(&b.inner.tensor);
    let a_slots: Vec<IndexSlot<'static>> = a_slots_data
        .iter()
        .map(|(label, pos)| {
            let static_label = unsafe { extend_str_lifetime(label) };
            IndexSlot::new(ctx.var(static_label), *pos)
        })
        .collect();
    let b_slots: Vec<IndexSlot<'static>> = b_slots_data
        .iter()
        .map(|(label, pos)| {
            let static_label = unsafe { extend_str_lifetime(label) };
            IndexSlot::new(ctx.var(static_label), *pos)
        })
        .collect();
    let a_rebuilt = Tensor::new(
        Symbol::new(unsafe { extend_str_lifetime(&a_name) }),
        a_slots,
    )
    .with_symmetry(a_sym);
    let b_rebuilt = Tensor::new(
        Symbol::new(unsafe { extend_str_lifetime(&b_name) }),
        b_slots,
    )
    .with_symmetry(b_sym);

    let result = contract(ctx, &a_rebuilt, &b_rebuilt);
    match result {
        Contracted::Product(p) => {
            let mut out: Vec<PyTensor> = Vec::with_capacity(p.factors.len());
            for factor in &p.factors {
                let (name, sym, slots) = snapshot(factor);
                out.push(rebuild_tensor(&name, sym, &slots)?);
            }
            let list = PyList::new(py, out)?;
            let tuple = ("product", list.into_any()).into_pyobject(py)?;
            Ok(tuple.into_any())
        }
        Contracted::Scalar(atom) => {
            let s = atom.to_string();
            let tuple = ("scalar", s).into_pyobject(py)?;
            Ok(tuple.into_any())
        }
    }
}

/// Return the symmetrisation sign of a tensor (+1 or -1).
///
/// For `symmetry="none"` and `"symmetric"` this is always +1. For
/// `"antisymmetric"` it returns the parity of the slot-sorting permutation.
#[pyfunction]
pub fn tensor_symmetrise_sign(tensor: &PyTensor) -> i64 {
    symmetrise_sign(&tensor.inner.tensor)
}
