//! Extension integration-rule families (0.27.1).
//!
//! Each family's specs live in their own `rules_ext_<family>.rs` sibling
//! module so per-family work lands in independent files. New family modules
//! register their specs in [`specs`] below and are compiled into the same
//! [`super::rules::IntegralRuleTable`] as the base library.

use super::rules::RuleSpec;

/// All extension specs, concatenated in family order.
pub(crate) fn specs() -> Vec<RuleSpec> {
    Vec::new()
}
