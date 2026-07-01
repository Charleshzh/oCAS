//! Function registry for user-defined external functions.
//!
//! The [`FunctionMap`] allows registering custom functions that can be
//! called during expression evaluation via [`ExternalFun`](crate::instruction::Instr::ExternalFun)
//! instructions.

use std::collections::HashMap;

use crate::domain::EvaluationDomain;

type ExternalFn<T> = Box<dyn Fn(&[T]) -> T + Send + Sync>;

/// A map of named functions that can be called during evaluation.
pub struct FunctionMap<T: EvaluationDomain> {
    entries: Vec<(String, FunctionEntry<T>)>,
    name_to_idx: HashMap<String, usize>,
    aliases: HashMap<String, String>,
}

/// A registered external function.
pub struct FunctionEntry<T: EvaluationDomain> {
    /// Number of arguments the function expects.
    pub arity: usize,
    func: ExternalFn<T>,
}

impl<T: EvaluationDomain> FunctionMap<T> {
    /// Create an empty function map.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            name_to_idx: HashMap::new(),
            aliases: HashMap::new(),
        }
    }

    /// Register a function with the given name and arity.
    pub fn register(&mut self, name: &str, arity: usize, func: ExternalFn<T>) {
        let idx = self.entries.len();
        self.entries
            .push((name.to_string(), FunctionEntry { arity, func }));
        self.name_to_idx.insert(name.to_string(), idx);
    }

    /// Register an alias for a function name.
    pub fn register_alias(&mut self, alias: &str, canonical: &str) {
        self.aliases
            .insert(alias.to_string(), canonical.to_string());
    }

    /// Look up a function by name (resolving aliases and case).
    pub fn resolve(&self, name: &str) -> Option<&FunctionEntry<T>> {
        self.resolve_idx(name).map(|idx| &self.entries[idx].1)
    }

    /// Get the index of a function by name.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.resolve_idx(name)
    }

    fn resolve_idx(&self, name: &str) -> Option<usize> {
        if let Some(idx) = self.name_to_idx.get(name) {
            return Some(*idx);
        }
        let lower = name.to_lowercase();
        if let Some(idx) = self.name_to_idx.get(&lower) {
            return Some(*idx);
        }
        if let Some(canonical) = self.aliases.get(name) {
            return self.name_to_idx.get(canonical.as_str()).copied();
        }
        if let Some(canonical) = self.aliases.get(&lower) {
            return self.name_to_idx.get(canonical.as_str()).copied();
        }
        None
    }

    /// Call a function by its index in the map.
    pub fn call_by_index(&self, idx: usize, args: &[T]) -> Option<T> {
        self.entries.get(idx).map(|(_, entry)| (entry.func)(args))
    }

    /// Return the number of registered functions.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Return true if no functions are registered.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

impl<T: EvaluationDomain> Default for FunctionMap<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn register_and_resolve() {
        let mut map = FunctionMap::<f64>::new();
        map.register("square", 1, Box::new(|args| args[0] * args[0]));
        assert!(map.resolve("square").is_some());
        assert!(map.resolve("Square").is_some());
        assert!(map.resolve("unknown").is_none());
    }

    #[test]
    fn alias_resolution() {
        let mut map = FunctionMap::<f64>::new();
        map.register("log", 1, Box::new(|args| args[0].ln()));
        map.register_alias("ln", "log");
        assert!(map.resolve("ln").is_some());
        assert!(map.resolve("Ln").is_some());
    }

    #[test]
    fn call_by_index() {
        let mut map = FunctionMap::<f64>::new();
        map.register("square", 1, Box::new(|args| args[0] * args[0]));
        let result = map.call_by_index(0, &[3.0]).unwrap();
        assert!((result - 9.0).abs() < 1e-10);
    }

    #[test]
    fn empty_map() {
        let map = FunctionMap::<f64>::new();
        assert!(map.is_empty());
        assert_eq!(map.len(), 0);
    }
}
