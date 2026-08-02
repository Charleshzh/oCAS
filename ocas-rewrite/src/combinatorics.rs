//! Combinatorial utilities: partition enumeration for `Transformer::Partition`.
//!
//! Implements the multiset partition algorithm used by Symbolica's
//! `Transformer::Partition`. Given a list of elements and a list of bins
//! (each with a name and capacity), enumerates all ways to distribute the
//! elements into the bins, with a multinomial coefficient prefactor.

use std::collections::HashMap;
use std::hash::Hash;

/// A single partition solution: a coefficient and an assignment of elements
/// to named bins.
#[derive(Debug, Clone)]
pub struct PartitionSolution<T: Clone + Ord + Hash, B: Clone + Eq + Hash> {
    /// Multinomial coefficient prefactor.
    pub coefficient: usize,
    /// Bin assignments (bin name → elements in that bin).
    pub bins: Vec<(B, Vec<T>)>,
}

/// Enumerate all ways to partition `elements` into `bins` (each is a
/// `(name, capacity)` pair).
///
/// * `fill_last` — if there are more elements than total bin capacity, the
///   surplus is absorbed into the last bin.
/// * `repeat` — if there are enough elements, the bin pattern is repeated
///   until all elements are consumed.
pub fn partitions<T: Clone + Ord + Hash, B: Clone + Ord + Hash>(
    elements: &[T],
    bins: &[(B, usize)],
    fill_last: bool,
    repeat: bool,
) -> Vec<PartitionSolution<T, B>> {
    if bins.is_empty() || elements.is_empty() {
        return Vec::new();
    }

    let bin_sum: usize = bins.iter().map(|b| b.1).sum();
    let total = elements.len();

    match total.cmp(&bin_sum) {
        std::cmp::Ordering::Less => return Vec::new(),
        std::cmp::Ordering::Equal => {}
        std::cmp::Ordering::Greater => {
            if !fill_last && (!repeat || !total.is_multiple_of(bin_sum)) {
                return Vec::new();
            }
        }
    }

    // Group equal elements.
    let mut element_groups: HashMap<T, usize> = HashMap::new();
    for e in elements {
        *element_groups.entry(e.clone()).or_insert(0) += 1;
    }
    let mut element_counts: Vec<(T, usize)> = element_groups.into_iter().collect();
    element_counts.sort_by(|a, b| a.0.cmp(&b.0));

    // Prepare bins.
    let mut sorted_bins = bins.to_vec();
    if fill_last {
        let last = sorted_bins.last_mut().unwrap();
        last.1 += total - bin_sum;
    }
    if repeat {
        for _ in 1..(total / bin_sum) {
            sorted_bins.extend_from_slice(bins);
        }
    }
    // Sort largest capacity first.
    sorted_bins.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let mut results: Vec<PartitionSolution<T, B>> = Vec::new();
    fill_rec(
        &sorted_bins,
        &mut element_counts,
        &mut Vec::new(),
        &mut Vec::new(),
        &mut Vec::new(),
        &mut results,
    );

    // Compute prefactors.
    for sol in &mut results {
        let mut coeff = 1usize;
        let mut counter = Vec::new();

        // Multinomial: for each element type, multinomial of its distribution across bins.
        for (elem, _total) in &element_counts {
            counter.clear();
            for (_, bin) in &sol.bins {
                let c = bin.iter().filter(|be| *be == elem).count();
                if c > 0 {
                    counter.push(c);
                }
            }
            coeff *= multinomial(&counter);
        }

        // Divide by factorial of identical bin configurations.
        let mut bin_groups: HashMap<&(B, Vec<T>), usize> = HashMap::new();
        for named_bin in &sol.bins {
            *bin_groups.entry(named_bin).or_insert(0) += 1;
        }
        for (_, count) in bin_groups {
            coeff /= factorial(count);
        }

        sol.coefficient = coeff;
    }

    results
}

/// Fill bins recursively.
fn fill_rec<T: Clone + Ord + Hash, B: Clone + Eq + Hash>(
    bins: &[(B, usize)],
    elem_counts: &mut [(T, usize)],
    single_buf: &mut Vec<T>,
    single_results: &mut Vec<Vec<T>>,
    accum: &mut Vec<(B, Vec<T>)>,
    results: &mut Vec<PartitionSolution<T, B>>,
) {
    if bins.is_empty() {
        if elem_counts.iter().all(|(_, c)| *c == 0) {
            results.push(PartitionSolution {
                coefficient: 1,
                bins: accum.clone(),
            });
        }
        return;
    }

    let (bin_id, bin_len) = &bins[0];
    let bin_id = bin_id.clone();
    let bin_len = *bin_len;

    // Enumerate all ways to fill this bin.
    single_results.clear();
    fill_bin(bin_len, elem_counts, single_buf, single_results);

    for fill in std::mem::take(single_results) {
        // Descending order check for duplicate bin names.
        if let Some(last) = accum.last()
            && last.0 == bin_id
            && fill.len() == last.1.len()
            && fill < last.1
        {
            continue;
        }

        // Deduct element counts.
        for x in &fill {
            if let Some((_, c)) = elem_counts.iter_mut().find(|(e, _)| *e == *x) {
                *c -= 1;
            }
        }

        accum.push((bin_id.clone(), fill.clone()));
        fill_rec(
            &bins[1..],
            elem_counts,
            single_buf,
            single_results,
            accum,
            results,
        );
        accum.pop();

        // Restore counts.
        for x in &fill {
            if let Some((_, c)) = elem_counts.iter_mut().find(|(e, _)| *e == *x) {
                *c += 1;
            }
        }
    }
}

/// Fill a single bin of `len` slots from available element counts.
fn fill_bin<T: Clone>(
    len: usize,
    elem_counts: &mut [(T, usize)],
    accum: &mut Vec<T>,
    results: &mut Vec<Vec<T>>,
) {
    if len == 0 {
        results.push(accum.clone());
        return;
    }
    let n = elem_counts.len();
    for i in 0..n {
        let count = elem_counts[i].1;
        if count > 0 {
            elem_counts[i].1 = count - 1;
            let name = elem_counts[i].0.clone();
            accum.push(name);
            fill_bin(len - 1, &mut elem_counts[i..], accum, results);
            accum.pop();
            elem_counts[i].1 = count;
        }
    }
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

fn multinomial(counts: &[usize]) -> usize {
    let total: usize = counts.iter().sum();
    let mut result = factorial(total);
    for &c in counts {
        result /= factorial(c);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_basic_exact() {
        let elements = vec![1i32, 3, 2, 3, 1];
        let bins = vec![('f', 2), ('g', 2), ('f', 1)];
        let sols = partitions(&elements, &bins, false, false);
        // There should be multiple valid ways to partition.
        assert!(!sols.is_empty(), "should produce at least one partition");
        for sol in &sols {
            let total: usize = sol.bins.iter().map(|(_, v)| v.len()).sum();
            assert_eq!(total, 5, "each solution must partition all 5 elements");
            assert!(sol.coefficient > 0, "coefficient must be positive");
        }
    }

    #[test]
    fn empty_inputs() {
        assert!(partitions::<i32, char>(&[], &[('a', 1)], false, false).is_empty());
        assert!(partitions::<i32, char>(&[1], &[], false, false).is_empty());
    }
}
