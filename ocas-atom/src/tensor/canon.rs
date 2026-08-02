//! Tensor expression canonicalisation via graph isomorphism.
//!
//! Encodes a tensor-product expression (Mul of Fun nodes) into a graph whose
//! vertex colours represent tensor heads / index slots and whose edges
//! represent argument positions and contractions.  The graph-isomorphism
//! engine [`super::graph`] then computes a canonical labelling, and the
//! result is reconstructed as a normalised tensor expression with renamed
//! dummy indices and reordered symmetric slots.

use std::collections::HashMap;

use crate::{Atom, AtomArena, AtomNode, Symbol};

use super::graph::{CanonicalForm, Graph};
use super::spec::TensorRegistry;

/// Error while canonicalising a tensor expression.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TensorCanonError {
    ContractedMoreThanOnce(Symbol),
    BadContraction(Symbol),
    NotATensor(Symbol),
    InconsistentOpenIndices,
    UnsupportedPower,
}

/// Result of canonicalising a tensor expression.
#[derive(Debug, Clone)]
pub struct CanonicalTensor<'a> {
    pub canonical_form: Atom<'a>,
    pub external_indices: Vec<Atom<'a>>,
    pub dummy_indices: Vec<Atom<'a>>,
}

// =========================================================================
// Public API
// =========================================================================

/// Canonicalise a tensor expression.
pub fn canonicalize_tensors<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<CanonicalTensor<'a>, TensorCanonError> {
    eprintln!("[canon-top] expr={}", expr);
    match expr.node() {
        AtomNode::Add(terms) => {
            let mut canon_terms: Vec<Atom<'a>> = Vec::new();
            let mut first_external: Option<Vec<Atom<'a>>> = None;
            let mut all_dummies: Vec<Atom<'a>> = Vec::new();

            for term in terms.iter() {
                let ct = canonicalize_single_term(ctx, *term, registry)?;
                match &first_external {
                    None => first_external = Some(ct.external_indices.clone()),
                    Some(ext) if *ext != ct.external_indices => {
                        return Err(TensorCanonError::InconsistentOpenIndices);
                    }
                    _ => {}
                }
                all_dummies.extend(ct.dummy_indices);
                canon_terms.push(ct.canonical_form);
            }

            let canonical_form = if canon_terms.len() == 1 {
                canon_terms.pop().unwrap()
            } else {
                ctx.add(&canon_terms)
            };
            Ok(CanonicalTensor {
                canonical_form,
                external_indices: first_external.unwrap_or_default(),
                dummy_indices: all_dummies,
            })
        }
        _ => canonicalize_single_term(ctx, expr, registry),
    }
}

fn canonicalize_single_term<'a>(
    ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<CanonicalTensor<'a>, TensorCanonError> {
    // Fast path: single tensor with all-symmetric slots. Sort args
    // alphabetically so the result is input-order-independent.
    #[allow(clippy::collapsible_if)]
    if let AtomNode::Fun(name, args) = expr.node() {
        eprintln!("[canon-1] Fun match: name={}", name.as_str());
        if let Some(spec) = registry.spec(*name) {
            eprintln!(
                "[canon-2] spec found: sym_subsets={}",
                spec.symmetric_subsets.len()
            );
            let all_symmetric = !spec.symmetric_subsets.is_empty()
                && spec.antisymmetric_subsets.is_empty()
                && spec.symmetric_subsets[0].len() == args.len();
            eprintln!("[canon-3] all_symmetric={}", all_symmetric);
            if all_symmetric {
                eprintln!("[canon-fast]");
                let mut sorted: Vec<Atom<'a>> = args.to_vec();
                sorted.sort_by_key(|a| match a.node() {
                    AtomNode::Var(s) => s.as_str().to_string(),
                    _ => a.to_string(),
                });
                let result = ctx.fun(name.as_str(), &sorted);
                return Ok(CanonicalTensor {
                    canonical_form: result,
                    external_indices: sorted,
                    dummy_indices: Vec::new(),
                });
            }
        }
    }

    let (g, head_nodes, slot_labels) = tensor_to_graph(ctx, expr, registry)?;
    let cf = g.canonize();
    reconstruct(ctx, &cf, &head_nodes, &slot_labels, registry)
}

// =========================================================================
// Graph encoding
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum TgNode {
    Head(u64),
    Slot(u64),
    Scalar(u64),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum TgEdge {
    HeadToSlot(usize, u8),
    Contraction(u64),
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct HeadInfo {
    symbol: Symbol,
    slot_count: usize,
    head_v: usize,
    slot_verts: Vec<usize>,
}

/// slot_vertex_index → original label Atom for all slots.
#[allow(clippy::type_complexity)]
fn tensor_to_graph<'a>(
    _ctx: &'a AtomArena<'a>,
    expr: Atom<'a>,
    registry: &TensorRegistry,
) -> Result<
    (
        Graph<TgNode, usize, TgEdge>,
        Vec<HeadInfo>,
        HashMap<usize, Atom<'a>>,
    ),
    TensorCanonError,
> {
    let mut g: Graph<TgNode, usize, TgEdge> = Graph::new();
    let mut heads: Vec<HeadInfo> = Vec::new();
    let mut index_uses: HashMap<Atom<'a>, (Vec<usize>, usize)> = HashMap::new();
    let mut slot_labels: HashMap<usize, Atom<'a>> = HashMap::new();

    match expr.node() {
        AtomNode::Mul(factors) => {
            for f in factors.iter() {
                encode_factor(
                    *f,
                    registry,
                    &mut g,
                    &mut heads,
                    &mut index_uses,
                    &mut slot_labels,
                )?;
            }
        }
        _ => {
            encode_factor(
                expr,
                registry,
                &mut g,
                &mut heads,
                &mut index_uses,
                &mut slot_labels,
            )?;
        }
    }

    for (_label, (slot_verts, count)) in &index_uses {
        if *count > 2 {
            return Err(TensorCanonError::ContractedMoreThanOnce(Symbol::new(
                &_label.to_string(),
            )));
        }
        if *count == 2 {
            let group = registry.index_group(Symbol::new(&_label.to_string()));
            g.add_undirected_edge(slot_verts[0], slot_verts[1], TgEdge::Contraction(group));
        }
    }

    Ok((g, heads, slot_labels))
}

fn encode_factor<'a>(
    factor: Atom<'a>,
    registry: &TensorRegistry,
    g: &mut Graph<TgNode, usize, TgEdge>,
    heads: &mut Vec<HeadInfo>,
    index_uses: &mut HashMap<Atom<'a>, (Vec<usize>, usize)>,
    slot_labels: &mut HashMap<usize, Atom<'a>>,
) -> Result<(), TensorCanonError> {
    match factor.node() {
        AtomNode::Fun(name, args) => {
            let spec = registry
                .spec(*name)
                .ok_or(TensorCanonError::NotATensor(*name))?;
            let head_v = g.add_node(TgNode::Head(hash(name.as_str())), 0);
            let mut slot_verts = Vec::with_capacity(args.len());

            // Pre-sort symmetric slots by label so the graph encoding is
            // input-order-independent.  Non-symmetric slots keep their
            // original position.
            let mut sorted_args: Vec<(usize, Atom<'a>)> =
                args.iter().enumerate().map(|(i, a)| (i, *a)).collect();
            // Stable sort: symmetric slots sorted by label, others by pos.
            sorted_args.sort_by(|&(pa, aa), &(pb, ab)| {
                let ha = spec.is_slot_hidden(pa);
                let hb = spec.is_slot_hidden(pb);
                match (ha, hb) {
                    (true, true) => {
                        // Compare by Symbol name to avoid platform-dependent
                        // Display differences.
                        let sa = match aa.node() {
                            AtomNode::Var(s) => s.as_str(),
                            _ => "",
                        };
                        let sb = match ab.node() {
                            AtomNode::Var(s) => s.as_str(),
                            _ => "",
                        };
                        sa.cmp(sb)
                    }
                    (true, false) => std::cmp::Ordering::Less,
                    (false, true) => std::cmp::Ordering::Greater,
                    (false, false) => pa.cmp(&pb),
                }
            });

            for (sorted_idx, (orig_pos, arg)) in sorted_args.into_iter().enumerate() {
                let label = arg;
                let is_hidden = spec.is_slot_hidden(orig_pos);
                // Symmetric slots must share the same colour so the graph-
                // isomorphism engine can freely permute them, ensuring
                // canonicalise(g(a,b)) == canonicalise(g(b,a)).
                let slot_colour = if is_hidden {
                    TgNode::Slot(0)
                } else {
                    TgNode::Slot(hash(&label.to_string()))
                };
                let slot_v = g.add_node(slot_colour, 0);
                slot_labels.insert(slot_v, label);
                slot_verts.push(slot_v);

                // Use sorted index as edge pos so the graph encoding is
                // identical for symmetric inputs like g(a,b) and g(b,a).
                let edge_pos = if is_hidden { sorted_idx } else { orig_pos };
                let kind = if is_hidden {
                    TgEdge::HeadToSlot(edge_pos, 0)
                } else {
                    TgEdge::HeadToSlot(edge_pos, 1)
                };
                g.add_directed_edge(head_v, slot_v, kind);

                let entry = index_uses.entry(label).or_insert_with(|| (Vec::new(), 0));
                entry.0.push(slot_v);
                entry.1 += 1;
            }

            heads.push(HeadInfo {
                symbol: *name,
                slot_count: args.len(),
                head_v,
                slot_verts,
            });
        }
        AtomNode::Pow(_, _) => return Err(TensorCanonError::UnsupportedPower),
        _ => {
            let h = hash(&factor.to_string());
            g.add_node(TgNode::Scalar(h), 0);
        }
    }
    Ok(())
}

// =========================================================================
// Reconstruction
// =========================================================================

#[allow(clippy::type_complexity, clippy::needless_range_loop)]
fn reconstruct<'a>(
    ctx: &'a AtomArena<'a>,
    cf: &CanonicalForm<TgNode, usize, TgEdge>,
    heads: &[HeadInfo],
    slot_labels: &HashMap<usize, Atom<'a>>,
    _registry: &TensorRegistry,
) -> Result<CanonicalTensor<'a>, TensorCanonError> {
    let cg = &cf.graph;
    let n = cg.node_count();

    // Map canonical → original vertex (vertex_map[pos] = original_vertex).
    let orig_of = &cf.vertex_map;

    // Find head→slot edges and contraction pairs in canonical graph.
    let mut slot_contraction: HashMap<usize, (usize, u64)> = HashMap::new();
    for v in 0..n {
        for ev in cg.edges_of(v) {
            if !ev.is_directed
                && let TgEdge::Contraction(g) = ev.data
            {
                slot_contraction.insert(v, (ev.neighbour, g));
            }
        }
    }

    // Assign canonical dummy names to each contraction pair.
    let mut group_counters: HashMap<u64, usize> = HashMap::new();
    // Key: (min(cv1, cv2), max(cv1, cv2))
    let mut pair_labels: HashMap<(usize, usize), Atom<'a>> = HashMap::new();

    for v in 0..n {
        for ev in cg.edges_of(v) {
            if !ev.is_directed
                && let TgEdge::Contraction(g) = ev.data
            {
                let a = v.min(ev.neighbour);
                let b = v.max(ev.neighbour);
                pair_labels.entry((a, b)).or_insert_with(|| {
                    let cnt = group_counters.entry(g).or_insert(0);
                    let label = if g == 0 {
                        ctx.var(&format!("d{}", cnt))
                    } else {
                        ctx.var(&format!("d{}_{}", g, cnt))
                    };
                    *cnt += 1;
                    label
                });
            }
        }
    }

    // Collect canonical heads and their slots.
    let mut canon_heads: Vec<(usize, &HeadInfo)> = Vec::new();
    let mut orig_to_head: HashMap<usize, &HeadInfo> = HashMap::new();
    for h in heads {
        orig_to_head.insert(h.head_v, h);
    }
    for v in 0..n {
        if let TgNode::Head(_) = cg.node_data(v) {
            let orig = orig_of[v];
            if let Some(h) = orig_to_head.get(&orig) {
                canon_heads.push((v, *h));
            }
        }
    }
    canon_heads.sort_by_key(|(v, _)| *v);

    // Build factors.
    let mut factors: Vec<Atom<'a>> = Vec::new();
    let mut all_dummies: Vec<Atom<'a>> = Vec::new();
    let mut external_indices: Vec<Atom<'a>> = Vec::new();

    for (can_head, h) in &canon_heads {
        // Gather slot info.
        let mut slot_infos: Vec<SlotInfo> = Vec::new();
        for ev in cg.edges_of(*can_head) {
            if ev.is_directed
                && ev.is_outgoing
                && let TgEdge::HeadToSlot(pos, hidden_flag) = ev.data
            {
                let partner = slot_contraction.get(&ev.neighbour).copied();
                slot_infos.push(SlotInfo {
                    orig_pos: pos,
                    hidden: hidden_flag == 0,
                    canon_slot_v: ev.neighbour,
                    partner_v: partner.map(|(p, _)| p),
                });
            }
        }

        // Sort: hidden (symmetric) slots by label string (deterministic
        // regardless of input order), visible slots by original position.
        slot_infos.sort_by(|a, b| match (a.hidden, b.hidden) {
            (true, true) => {
                let la = orig_of[a.canon_slot_v];
                let lb = orig_of[b.canon_slot_v];
                let sa = slot_labels
                    .get(&la)
                    .map(|x| x.to_string())
                    .unwrap_or_default();
                let sb = slot_labels
                    .get(&lb)
                    .map(|x| x.to_string())
                    .unwrap_or_default();
                sa.cmp(&sb)
            }
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            (false, false) => a.orig_pos.cmp(&b.orig_pos),
        });

        let mut args: Vec<Atom<'a>> = Vec::new();
        for si in &slot_infos {
            if let Some(pv) = si.partner_v {
                let a = si.canon_slot_v.min(pv);
                let b = si.canon_slot_v.max(pv);
                if let Some(label) = pair_labels.get(&(a, b)) {
                    args.push(*label);
                    if !all_dummies.contains(label) {
                        all_dummies.push(*label);
                    }
                } else {
                    args.push(ctx.var("?"));
                }
            } else {
                // External index: preserve original label from the graph encoding.
                let orig_slot = orig_of[si.canon_slot_v];
                if let Some(&orig_label) = slot_labels.get(&orig_slot) {
                    args.push(orig_label);
                    if !external_indices.contains(&orig_label) {
                        external_indices.push(orig_label);
                    }
                } else {
                    // Fallback: synthetic name.
                    let label = ctx.var(&format!("ext{}", external_indices.len()));
                    args.push(label);
                    if !external_indices.contains(&label) {
                        external_indices.push(label);
                    }
                }
            }
        }

        factors.push(ctx.fun(h.symbol.as_str(), &args));
    }

    let canonical_form = if factors.is_empty() {
        ctx.num(1)
    } else if factors.len() == 1 {
        factors.pop().unwrap()
    } else {
        ctx.mul(&factors)
    };

    Ok(CanonicalTensor {
        canonical_form,
        external_indices,
        dummy_indices: all_dummies,
    })
}

struct SlotInfo {
    orig_pos: usize,
    hidden: bool,
    canon_slot_v: usize,
    partner_v: Option<usize>,
}

fn hash(s: &str) -> u64 {
    use std::hash::Hasher;
    let mut h = std::collections::hash_map::DefaultHasher::new();
    std::hash::Hash::hash(&s, &mut h);
    h.finish()
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AtomArena;
    use crate::Symbol;
    use crate::tensor::spec::SymmetrySpec;
    use ocas_core::arena::Arena;

    #[test]
    fn canon_single_tensor_no_symmetry() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let mut reg = TensorRegistry::new();
        reg.register(Symbol::new("T"), SymmetrySpec::none());

        let i = ctx.var("i");
        let j = ctx.var("j");
        let t = ctx.fun("T", &[i, j]);
        let ct = canonicalize_tensors(&ctx, t, &reg).unwrap();
        let s = ct.canonical_form.to_string();
        assert!(s.contains("T"), "result: {s}");
    }

    #[test]
    fn canon_product_with_contraction() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let mut reg = TensorRegistry::new();
        reg.register(Symbol::new("T"), SymmetrySpec::none());
        reg.register(Symbol::new("U"), SymmetrySpec::none());

        let i = ctx.var("i");
        let j = ctx.var("j");
        let k = ctx.var("k");
        let t = ctx.fun("T", &[i, j]);
        let u = ctx.fun("U", &[j, k]);
        let prod = ctx.mul(&[t, u]);
        let ct = canonicalize_tensors(&ctx, prod, &reg).unwrap();
        let s = ct.canonical_form.to_string();
        // Should have at least one dummy and two tensors.
        assert!(s.contains("d0"), "expected dummy d0, got: {s}");
        assert!(s.contains('T') && s.contains('U'), "got: {s}");
        assert_eq!(ct.dummy_indices.len(), 1);
    }

    #[test]
    fn canon_symmetric_tensor_consistency() {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let mut reg = TensorRegistry::new();
        reg.register(Symbol::new("g"), SymmetrySpec::fully_symmetric(2));

        let a = ctx.var("a");
        let b = ctx.var("b");
        let g_ab = ctx.fun("g", &[a, b]);
        let g_ba = ctx.fun("g", &[b, a]);
        let ct1 = canonicalize_tensors(&ctx, g_ab, &reg).unwrap();
        let ct2 = canonicalize_tensors(&ctx, g_ba, &reg).unwrap();
        // Both should canonicalise to the same form.
        assert_eq!(
            ct1.canonical_form.to_string(),
            ct2.canonical_form.to_string(),
            "symmetric slots should canonicalise consistently"
        );
    }
}
