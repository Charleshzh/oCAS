//! Optional GPL-licensed backends for oCAS.
//!
//! This crate is a placeholder: it isolates any future GPL-3.0-or-later
//! dependencies (e.g. NTL, Singular, PARI/GP interfaces) from the default
//! LGPL-3.0-or-later build. It currently exposes no API and no backend.
//!
//! The top-level `ocas` facade intentionally does not expose a `gpl` feature:
//! enabling a GPL backend means depending on `ocas-gpl` directly once it has
//! content, so that the combined work is explicitly opt-in.
