//! Streaming evaluation over large datasets.
//!
//! [`StreamingEvaluator`] wraps an [`ExpressionEvaluator`] with reusable
//! buffers, so processing a stream of input rows uses constant memory
//! regardless of stream length. This mirrors the semantics of Symbolica's
//! `streaming.rs`: rows flow through the evaluator one at a time and
//! results are consumed by a sink callback.

use crate::domain::{EvaluationDomain, PowfExtension};
use crate::error::Result;
use crate::evaluator::ExpressionEvaluator;

/// Streaming evaluator with reusable internal buffers.
///
/// Created from a borrowed [`ExpressionEvaluator`]; all scratch memory
/// (parameter staging, evaluation stack, result buffer) is allocated
/// once up front and reused for every row.
///
/// # Example
///
/// ```ignore
/// let eval: ExpressionEvaluator<f64> = ExpressionEvaluator::compile(atom)?;
/// let mut stream = StreamingEvaluator::new(&eval);
/// let rows = (0..1_000_000).map(|i| [i as f64]);
/// let n = stream.for_each(rows, |results| {
///     // consume results (same buffer, valid only within the callback)
/// })?;
/// assert_eq!(n, 1_000_000);
/// ```
pub struct StreamingEvaluator<'a, T: EvaluationDomain> {
    evaluator: &'a ExpressionEvaluator<T>,
    params: Vec<T>,
    stack: Vec<T>,
    results: Vec<T>,
}

impl<T: EvaluationDomain + PowfExtension> StreamingEvaluator<'_, T> {
    /// Create a streaming evaluator, pre-allocating all buffers.
    pub fn new(evaluator: &ExpressionEvaluator<T>) -> StreamingEvaluator<'_, T> {
        StreamingEvaluator {
            evaluator,
            params: Vec::with_capacity(evaluator.param_count()),
            stack: Vec::with_capacity(evaluator.stack_size()),
            results: Vec::with_capacity(evaluator.result_count()),
        }
    }

    /// Process a stream of input rows, invoking `sink` with the result
    /// slice for each row.
    ///
    /// Each row must contain exactly
    /// [`param_count`](ExpressionEvaluator::param_count) values; the
    /// slice passed to `sink` contains
    /// [`result_count`](ExpressionEvaluator::result_count) values and is
    /// only valid for the duration of the callback. Returns the number
    /// of rows processed.
    ///
    /// Memory usage is constant: no allocation grows with the number of
    /// rows.
    ///
    /// # Errors
    ///
    /// Returns [`crate::EvaluationError`] if a row has the wrong arity
    /// or an arithmetic error occurs during evaluation.
    pub fn for_each<I, S, F>(&mut self, rows: I, mut sink: F) -> Result<usize>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[T]>,
        F: FnMut(&[T]),
    {
        let mut count = 0usize;
        for row in rows {
            let row = row.as_ref();
            self.params.clear();
            self.params.extend(row.iter().cloned());
            self.evaluator
                .evaluate_with_stack(&self.params, &mut self.stack, &mut self.results)?;
            sink(&self.results);
            count += 1;
        }
        Ok(count)
    }

    /// Process a chunk of rows and collect all results.
    ///
    /// Convenience method for bounded batches; prefer
    /// [`for_each`](StreamingEvaluator::for_each) for unbounded streams.
    /// Returns one result vector per row.
    ///
    /// # Errors
    ///
    /// Returns [`crate::EvaluationError`] if a row has the wrong arity
    /// or an arithmetic error occurs during evaluation.
    pub fn evaluate_chunk<S: AsRef<[T]>>(&mut self, rows: &[S]) -> Result<Vec<Vec<T>>> {
        let mut out = Vec::with_capacity(rows.len());
        self.for_each(rows, |results| out.push(results.to_vec()))?;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ocas_atom::AtomArena;
    use ocas_core::arena::Arena;

    fn build_eval() -> (Arena, ExpressionEvaluator<f64>) {
        let arena = Arena::new();
        let ctx = AtomArena::new(&arena);
        let sum = ctx.add(&[ctx.var("x"), ctx.var("y")]);
        let prod = ctx.mul(&[ctx.var("x"), ctx.var("y")]);
        let eval = ExpressionEvaluator::compile_multi(&[sum, prod]).unwrap();
        (arena, eval)
    }

    #[test]
    fn streaming_multi_output() {
        let (_arena, eval) = build_eval();
        let mut stream = StreamingEvaluator::new(&eval);
        let rows: Vec<[f64; 2]> = (0..100).map(|i| [i as f64, 2.0]).collect();
        let mut seen = Vec::new();
        let n = stream
            .for_each(&rows, |results| seen.push((results[0], results[1])))
            .unwrap();
        assert_eq!(n, 100);
        for (i, &(sum, prod)) in seen.iter().enumerate() {
            assert!((sum - (i as f64 + 2.0)).abs() < 1e-10);
            assert!((prod - (i as f64 * 2.0)).abs() < 1e-10);
        }
    }

    #[test]
    fn streaming_constant_memory_million_rows() {
        let (_arena, eval) = build_eval();
        let mut stream = StreamingEvaluator::new(&eval);

        // Warm up so buffers reach their steady-state capacity.
        let warm: Vec<[f64; 2]> = vec![[1.0, 2.0]; 10];
        stream.for_each(&warm, |_| {}).unwrap();
        let stack_cap = stream.stack.capacity();
        let results_cap = stream.results.capacity();
        let params_cap = stream.params.capacity();

        // One million rows produced lazily — no dataset allocation.
        let rows = (0..1_000_000u64).map(|i| [i as f64 % 100.0, 3.0]);
        let mut count = 0usize;
        let mut checksum = 0.0f64;
        let n = stream
            .for_each(rows, |results| {
                count += 1;
                checksum += results[0];
            })
            .unwrap();
        assert_eq!(n, 1_000_000);
        assert_eq!(count, 1_000_000);
        assert!(checksum > 0.0);

        // Buffer capacities unchanged: memory is constant in stream length.
        assert_eq!(stream.stack.capacity(), stack_cap);
        assert_eq!(stream.results.capacity(), results_cap);
        assert_eq!(stream.params.capacity(), params_cap);
    }

    #[test]
    fn streaming_wrong_arity_errors() {
        let (_arena, eval) = build_eval();
        let mut stream = StreamingEvaluator::new(&eval);
        let rows: Vec<Vec<f64>> = vec![vec![1.0], vec![2.0, 3.0]];
        assert!(stream.for_each(&rows, |_| {}).is_err());
    }

    #[test]
    fn streaming_evaluate_chunk() {
        let (_arena, eval) = build_eval();
        let mut stream = StreamingEvaluator::new(&eval);
        let rows: Vec<[f64; 2]> = vec![[1.0, 2.0], [3.0, 4.0]];
        let out = stream.evaluate_chunk(&rows).unwrap();
        assert_eq!(out.len(), 2);
        assert!((out[0][0] - 3.0).abs() < 1e-10);
        assert!((out[0][1] - 2.0).abs() < 1e-10);
        assert!((out[1][0] - 7.0).abs() < 1e-10);
        assert!((out[1][1] - 12.0).abs() < 1e-10);
    }
}
