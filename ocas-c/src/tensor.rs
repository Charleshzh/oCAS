//! C/C++ bindings for the basic tensor algebra module.
//!
//! Exposes [`ocas_atom::tensor`] — slot-based tensors with index variance and
//! symmetry, explicit contraction, and a symmetrisation sign — as opaque C
//! handles.
//!
//! # Slot string format
//!
//! Slots are passed as a single semicolon-separated string, each entry being
//! `label,position` where `position` is `upper` or `lower` (aliases: `up`,
//! `down`, `contravariant`, `covariant`). For example:
//!
//! - `"i,upper;j,lower"` — two slots `i^j`.
//! - `"a,lower;b,lower"` — two lower slots.
//!
//! # Contraction result convention
//!
//! [`ocas_tensor_contract`] writes the result into an [`OcasTensorContraction`]
//! struct. The `kind` field is `0` for `"product"` (free slots remain;
//! `tensors` holds the resulting tensor handles) or `1` for `"scalar"` (fully
//! contracted; `scalar_str` holds the string form of the atom). The caller
//! must free each tensor in `tensors` via [`ocas_tensor_free`] and release
//! `scalar_str` via [`ocas_string_free`].

#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::{CStr, CString, c_char, c_int};
use std::ptr;

use ocas_atom::tensor::{
    Contracted, IndexPosition, IndexSlot, Symmetry, Tensor, contract, symmetrise_sign,
};
use ocas_atom::{AtomArena, Symbol};
use ocas_core::arena::Arena;

use crate::error::{
    OCAS_ERROR_INVALID_ARGUMENT, OCAS_ERROR_NULL_POINTER, OCAS_ERROR_RUNTIME, OCAS_OK, set,
};

/// Opaque handle for a tensor.
#[repr(C)]
pub struct OcasTensor {
    _private: [u8; 0],
}

/// A result of [`ocas_tensor_contract`].
#[repr(C)]
pub struct OcasTensorContraction {
    /// `0` = product (free slots remain), `1` = scalar (fully contracted).
    pub kind: c_int,
    /// Pointer to an array of `OcasTensor*` handles (valid when `kind == 0`).
    /// May be `NULL` if `n_tensors == 0`.
    pub tensors: *mut *mut OcasTensor,
    /// Number of tensor handles in `tensors`.
    pub n_tensors: usize,
    /// Heap-allocated string form of the scalar atom (valid when `kind == 1`).
    /// `NULL` otherwise. Release with [`ocas_string_free`].
    pub scalar_str: *mut c_char,
}

// ------------------------------------------------------------------
//  Arena management
// ------------------------------------------------------------------

/// Extend a string's lifetime to `'static`. Safe because atoms never retain
/// borrows of the input string.
unsafe fn extend_str_lifetime(s: &str) -> &'static str {
    unsafe { std::mem::transmute::<&str, &'static str>(s) }
}

/// Storage behind an [`OcasTensor`]: a leaked arena pair recovered on drop.
struct TensorInner {
    arena_ptr: *mut Arena,
    ctx_ptr: *mut AtomArena<'static>,
    tensor: Tensor<'static>,
}

// SAFETY: matches ocas-c::expression::ExprBox.
unsafe impl Send for TensorInner {}

impl Drop for TensorInner {
    fn drop(&mut self) {
        // SAFETY: both pointers came from `Box::into_raw`; drop ctx first.
        unsafe {
            let _ = Box::from_raw(self.ctx_ptr);
            let _ = Box::from_raw(self.arena_ptr);
        }
    }
}

fn leak_arena_and_ctx() -> (*mut Arena, *mut AtomArena<'static>) {
    let arena_box: Box<Arena> = Box::new(Arena::new());
    let arena_ptr = Box::into_raw(arena_box);
    // SAFETY: arena_ptr outlives TensorInner; recovered in Drop.
    let arena_ref: &'static Arena = unsafe { &*arena_ptr };
    let ctx = AtomArena::new(arena_ref);
    let ctx_ptr = Box::into_raw(Box::new(ctx));
    (arena_ptr, ctx_ptr)
}

// ------------------------------------------------------------------
//  Parsing helpers
// ------------------------------------------------------------------

fn parse_position(s: &str) -> Option<IndexPosition> {
    match s.trim().to_ascii_lowercase().as_str() {
        "upper" | "up" | "contravariant" => Some(IndexPosition::Upper),
        "lower" | "down" | "covariant" => Some(IndexPosition::Lower),
        _ => None,
    }
}

fn parse_symmetry(s: &str) -> Option<Symmetry> {
    match s.trim().to_ascii_lowercase().as_str() {
        "none" | "" => Some(Symmetry::None),
        "symmetric" | "sym" => Some(Symmetry::Symmetric),
        "antisymmetric" | "antisym" | "skew" => Some(Symmetry::Antisymmetric),
        _ => None,
    }
}

/// Parse the slot string `"i,upper;j,lower"` into a list of (label, position).
fn parse_slots(slots_str: &str) -> Result<Vec<(String, IndexPosition)>, String> {
    let mut out = Vec::new();
    for entry in slots_str.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let mut parts = entry.split(',');
        let label = parts
            .next()
            .ok_or_else(|| format!("missing label in slot entry {entry:?}"))?
            .trim()
            .to_string();
        if label.is_empty() {
            return Err(format!("empty label in slot entry {entry:?}"));
        }
        let pos_str = parts
            .next()
            .ok_or_else(|| format!("missing position in slot entry {entry:?}"))?;
        let pos = parse_position(pos_str)
            .ok_or_else(|| format!("invalid position {pos_str:?} in slot entry {entry:?}"))?;
        if parts.next().is_some() {
            return Err(format!("extra fields in slot entry {entry:?}"));
        }
        out.push((label, pos));
    }
    Ok(out)
}

// ------------------------------------------------------------------
//  Opaque-handle helpers
// ------------------------------------------------------------------

fn tensor_ptr(t: Box<TensorInner>) -> *mut OcasTensor {
    Box::into_raw(t) as *mut OcasTensor
}

fn tensor_ref<'a>(t: *const OcasTensor) -> Option<&'a TensorInner> {
    if t.is_null() {
        return None;
    }
    Some(unsafe { &*(t as *const TensorInner) })
}

fn build_tensor_inner(
    name: &str,
    slots: &[(String, IndexPosition)],
    symmetry: Symmetry,
) -> Box<TensorInner> {
    let (arena_ptr, ctx_ptr) = leak_arena_and_ctx();
    let ctx: &'static AtomArena<'static> = unsafe { &*ctx_ptr };
    let static_name = unsafe { extend_str_lifetime(name) };
    let symbol = Symbol::new(static_name);
    let slots: Vec<IndexSlot<'static>> = slots
        .iter()
        .map(|(label, pos)| {
            let static_label = unsafe { extend_str_lifetime(label) };
            IndexSlot::new(ctx.var(static_label), *pos)
        })
        .collect();
    let tensor = Tensor::new(symbol, slots).with_symmetry(symmetry);
    Box::new(TensorInner {
        arena_ptr,
        ctx_ptr,
        tensor,
    })
}

/// Snapshot a tensor's name, symmetry, and slots as plain data.
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

// ------------------------------------------------------------------
//  C API
// ------------------------------------------------------------------

/// Create a tensor from a name and a slots string (see the
/// [module docs](self) for the format). The optional symmetry is one of
/// `"none"`, `"symmetric"`, `"antisymmetric"`; pass `NULL` for `"none"`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_create(
    name: *const c_char,
    slots: *const c_char,
    symmetry: *const c_char,
    err: *mut c_int,
) -> *mut OcasTensor {
    crate::error::clear();
    if name.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "tensor name string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    if slots.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "slots string is null");
        crate::error::write_last_code(err);
        return ptr::null_mut();
    }
    let name_str = unsafe { CStr::from_ptr(name) };
    let name_str = match name_str.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(
                OCAS_ERROR_INVALID_ARGUMENT,
                "name string is not valid UTF-8",
            );
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let slots_str = unsafe { CStr::from_ptr(slots) };
    let slots_str = match slots_str.to_str() {
        Ok(s) => s,
        Err(_) => {
            set(
                OCAS_ERROR_INVALID_ARGUMENT,
                "slots string is not valid UTF-8",
            );
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let sym = if symmetry.is_null() {
        Symmetry::None
    } else {
        let s = unsafe { CStr::from_ptr(symmetry) };
        match s.to_str() {
            Ok(ss) => match parse_symmetry(ss) {
                Some(s) => s,
                None => {
                    set(
                        OCAS_ERROR_INVALID_ARGUMENT,
                        "invalid symmetry string (use none/symmetric/antisymmetric)",
                    );
                    crate::error::write_last_code(err);
                    return ptr::null_mut();
                }
            },
            Err(_) => {
                set(
                    OCAS_ERROR_INVALID_ARGUMENT,
                    "symmetry string is not valid UTF-8",
                );
                crate::error::write_last_code(err);
                return ptr::null_mut();
            }
        }
    };
    let parsed = match parse_slots(slots_str) {
        Ok(v) => v,
        Err(msg) => {
            set(OCAS_ERROR_INVALID_ARGUMENT, &msg);
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let inner = build_tensor_inner(name_str, &parsed, sym);
    crate::error::write_last_code(err);
    tensor_ptr(inner)
}

/// Free a tensor handle. Safe to call with `NULL`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_free(t: *mut OcasTensor) {
    if t.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(t as *mut TensorInner));
    }
}

/// Return the tensor's name as a heap-allocated string. The caller must
/// release it with [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_name(t: *const OcasTensor, err: *mut c_int) -> *mut c_char {
    crate::error::clear();
    let inner = match tensor_ref(t) {
        Some(i) => i,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "tensor handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    match CString::new(inner.tensor.name().as_str()) {
        Ok(cs) => {
            crate::error::write_last_code(err);
            cs.into_raw()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "tensor name contains a NUL byte");
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Return the tensor's arity (number of slots), or `0` on a null handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_rank(t: *const OcasTensor) -> usize {
    match tensor_ref(t) {
        Some(i) => i.tensor.rank(),
        None => 0,
    }
}

/// Return the tensor's symmetry as an integer code, or `-1` on a null handle:
/// `0` = none, `1` = symmetric, `2` = antisymmetric.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_symmetry(t: *const OcasTensor) -> c_int {
    match tensor_ref(t) {
        Some(i) => match i.tensor.symmetry() {
            Symmetry::None => 0,
            Symmetry::Symmetric => 1,
            Symmetry::Antisymmetric => 2,
        },
        None => -1,
    }
}

/// Return a heap-allocated string representation of the tensor rendered as
/// an atom `name(slot, slot, ...)`. The caller must release it with
/// [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_to_string(t: *const OcasTensor, err: *mut c_int) -> *mut c_char {
    crate::error::clear();
    let inner = match tensor_ref(t) {
        Some(i) => i,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "tensor handle is null");
            crate::error::write_last_code(err);
            return ptr::null_mut();
        }
    };
    let ctx: &'static AtomArena<'static> = unsafe { &*inner.ctx_ptr };
    let s = inner.tensor.to_atom(ctx).to_string();
    match CString::new(s) {
        Ok(cs) => {
            crate::error::write_last_code(err);
            cs.into_raw()
        }
        Err(_) => {
            set(OCAS_ERROR_RUNTIME, "tensor string contains a NUL byte");
            crate::error::write_last_code(err);
            ptr::null_mut()
        }
    }
}

/// Return the symmetrisation sign of the tensor (+1 or -1), or `0` on a null
/// handle.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_symmetrise_sign(t: *const OcasTensor) -> i64 {
    match tensor_ref(t) {
        Some(i) => symmetrise_sign(&i.tensor),
        None => 0,
    }
}

/// Contract two tensors by summing over shared dummy indices.
///
/// On success `out` is filled. When `out.kind == 0` (product), `out.tensors`
/// holds `out.n_tensors` independent tensor handles that the caller must free
/// via [`ocas_tensor_free`]. When `out.kind == 1` (scalar), `out.scalar_str`
/// holds a heap-allocated string that must be released via
/// [`ocas_string_free`].
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_contract(
    a: *const OcasTensor,
    b: *const OcasTensor,
    out: *mut OcasTensorContraction,
    err: *mut c_int,
) -> c_int {
    crate::error::clear();
    if out.is_null() {
        set(OCAS_ERROR_NULL_POINTER, "output pointer is null");
        crate::error::write_last_code(err);
        return OCAS_ERROR_NULL_POINTER;
    }
    let a_inner = match tensor_ref(a) {
        Some(i) => i,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "first tensor handle is null");
            crate::error::write_last_code(err);
            return OCAS_ERROR_NULL_POINTER;
        }
    };
    let b_inner = match tensor_ref(b) {
        Some(i) => i,
        None => {
            set(OCAS_ERROR_NULL_POINTER, "second tensor handle is null");
            crate::error::write_last_code(err);
            return OCAS_ERROR_NULL_POINTER;
        }
    };

    // Allocate a shared arena for the contraction computation; dropped before
    // returning. The result tensors are rebuilt into independent arenas.
    let (arena_ptr, ctx_ptr) = leak_arena_and_ctx();
    struct DropGuard(*mut Arena, *mut AtomArena<'static>);
    impl Drop for DropGuard {
        fn drop(&mut self) {
            unsafe {
                let _ = Box::from_raw(self.1);
                let _ = Box::from_raw(self.0);
            }
        }
    }
    let _guard = DropGuard(arena_ptr, ctx_ptr);
    let ctx: &'static AtomArena<'static> = unsafe { &*ctx_ptr };

    // Rebuild a and b into the shared arena.
    let (a_name, a_sym, a_slots_data) = snapshot(&a_inner.tensor);
    let (b_name, b_sym, b_slots_data) = snapshot(&b_inner.tensor);
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
            let mut handles: Vec<*mut OcasTensor> = Vec::with_capacity(p.factors.len());
            for factor in &p.factors {
                let (name, sym, slots) = snapshot(factor);
                let inner = build_tensor_inner(&name, &slots, sym);
                handles.push(tensor_ptr(inner));
            }
            let n = handles.len();
            // `into_boxed_slice` guarantees len == capacity, so the later
            // `Vec::from_raw_parts(ptr, n, n)` recovery is sound.
            let mut boxed = handles.into_boxed_slice();
            let raw = boxed.as_mut_ptr();
            std::mem::forget(boxed);
            unsafe {
                ptr::write(
                    out,
                    OcasTensorContraction {
                        kind: 0,
                        tensors: raw,
                        n_tensors: n,
                        scalar_str: ptr::null_mut(),
                    },
                );
            }
            crate::error::write_last_code(err);
            crate::error::OCAS_OK
        }
        Contracted::Scalar(atom) => {
            let s = atom.to_string();
            match CString::new(s) {
                Ok(cs) => {
                    let raw = cs.into_raw();
                    unsafe {
                        ptr::write(
                            out,
                            OcasTensorContraction {
                                kind: 1,
                                tensors: ptr::null_mut(),
                                n_tensors: 0,
                                scalar_str: raw,
                            },
                        );
                    }
                    crate::error::write_last_code(err);
                    OCAS_OK
                }
                Err(_) => {
                    set(OCAS_ERROR_RUNTIME, "scalar string contains a NUL byte");
                    crate::error::write_last_code(err);
                    OCAS_ERROR_RUNTIME
                }
            }
        }
    }
}

/// Free the `tensors` array (but NOT the individual tensor handles) returned
/// by [`ocas_tensor_contract`] when `kind == 0`. Each tensor handle must be
/// freed separately via [`ocas_tensor_free`]. Safe to call with `NULL`.
#[unsafe(no_mangle)]
pub extern "C" fn ocas_tensor_contraction_free(c: *mut OcasTensorContraction) {
    if c.is_null() {
        return;
    }
    unsafe {
        let c_ref = &mut *c;
        if !c_ref.tensors.is_null() && c_ref.n_tensors > 0 {
            let _ = Vec::from_raw_parts(c_ref.tensors, c_ref.n_tensors, c_ref.n_tensors);
        }
        c_ref.tensors = ptr::null_mut();
        c_ref.n_tensors = 0;
        if !c_ref.scalar_str.is_null() {
            crate::expression::ocas_string_free(c_ref.scalar_str);
            c_ref.scalar_str = ptr::null_mut();
        }
    }
}
