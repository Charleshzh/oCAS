//! Arena allocator for expression nodes.
//!
//! oCAS uses a bump allocator to store expression sub-nodes. This avoids the
//! per-node allocation overhead of `Box` or `Rc` and improves cache locality.
//! When the arena is dropped, the entire tree is freed at once.
//!
//! # Current limitations
//!
//! The 0.1.0 `Arena` does **not** run destructors for allocated values. It is
//! therefore only safe to store `Copy` types or types that do not own resources
//! requiring explicit cleanup. This restriction will be lifted once expression
//! trees need to store owned strings or other `Drop` types.

use std::alloc::{self, Layout};
use std::cell::RefCell;
use std::marker::PhantomData;
use std::mem;
use std::ptr::NonNull;

/// Default block size for new arena chunks.
const DEFAULT_BLOCK_SIZE: usize = 64 * 1024;

/// A bump-allocated region of memory.
///
/// Values allocated in an `Arena` are tied to its lifetime and must not
/// outlive it. The public API enforces this with borrow checker lifetimes.
///
/// # Type safety note
///
/// `Arena` does not run destructors. Only store `Copy` or otherwise
/// non-owning values until drop support is added.
pub struct Arena {
    chunks: RefCell<Vec<Chunk>>,
    block_size: usize,
}

impl Arena {
    /// Create a new arena with the default block size.
    pub fn new() -> Self {
        Self {
            chunks: RefCell::new(Vec::new()),
            block_size: DEFAULT_BLOCK_SIZE,
        }
    }

    /// Create a new arena with a custom initial block size.
    pub fn with_capacity(block_size: usize) -> Self {
        Self {
            chunks: RefCell::new(Vec::new()),
            block_size,
        }
    }

    /// Allocate a value in the arena, constructing it inside `init`, and return
    /// a mutable reference tied to `self`.
    ///
    /// The closure form avoids any ambiguity about when mutation of the arena
    /// occurs. The returned `&mut T` is unique because `alloc_raw` advances the
    /// arena offset for each allocation via interior mutability.
    ///
    /// # Panics
    ///
    /// Panics if the requested layout has size zero.
    #[allow(clippy::mut_from_ref)]
    pub fn allocate_with<T>(&self, init: impl FnOnce() -> T) -> &mut T {
        let layout = Layout::new::<T>();
        assert!(
            layout.size() > 0,
            "cannot allocate zero-sized types in Arena"
        );
        let ptr = self.alloc_raw(layout);

        // SAFETY: `ptr` is non-null and properly aligned for `T`.
        unsafe {
            let typed = ptr.as_ptr().cast::<T>();
            typed.write(init());
            &mut *typed
        }
    }

    fn alloc_raw(&self, layout: Layout) -> NonNull<u8> {
        let mut chunks = self.chunks.borrow_mut();

        // Try to allocate from the current chunk.
        if let Some(chunk) = chunks.last_mut()
            && let Some(ptr) = chunk.try_alloc(layout)
        {
            return ptr;
        }

        // Need a new chunk. Use at least the requested size.
        let size = layout.size().max(self.block_size);
        let mut new_chunk = Chunk::new(size);
        let ptr = new_chunk
            .try_alloc(layout)
            .expect("new chunk should fit any layout up to its size");
        chunks.push(new_chunk);
        ptr
    }
}

impl Default for Arena {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for Arena {
    fn drop(&mut self) {
        // 0.1.0: destructors for allocated values are intentionally not called.
        // Only `Copy`/non-owning values may be stored.
    }
}

struct Chunk {
    memory: NonNull<u8>,
    size: usize,
    offset: usize,
}

impl Chunk {
    fn new(size: usize) -> Self {
        let layout =
            Layout::from_size_align(size, mem::align_of::<usize>()).expect("invalid chunk layout");
        // SAFETY: layout is non-zero and properly aligned.
        let memory = unsafe { NonNull::new_unchecked(alloc::alloc(layout)) };
        Self {
            memory,
            size,
            offset: 0,
        }
    }

    fn try_alloc(&mut self, layout: Layout) -> Option<NonNull<u8>> {
        let aligned_offset = align_up(self.offset, layout.align());
        let end = aligned_offset.checked_add(layout.size())?;
        if end > self.size {
            return None;
        }

        // SAFETY: offset is within the allocated block and aligned.
        let ptr = unsafe { NonNull::new_unchecked(self.memory.as_ptr().add(aligned_offset)) };
        self.offset = end;
        Some(ptr)
    }
}

impl Drop for Chunk {
    fn drop(&mut self) {
        let layout = Layout::from_size_align(self.size, mem::align_of::<usize>())
            .expect("invalid chunk layout");
        // SAFETY: `memory` was allocated with this layout.
        unsafe {
            alloc::dealloc(self.memory.as_ptr(), layout);
        }
    }
}

fn align_up(offset: usize, align: usize) -> usize {
    assert!(align.is_power_of_two(), "alignment must be a power of two");
    (offset + align - 1) & !(align - 1)
}

/// An owned expression that keeps its arena alive.
pub struct OwnedExpr<T> {
    #[allow(dead_code)]
    arena: Box<Arena>,
    root: *mut T,
    _marker: PhantomData<T>,
}

impl<T> OwnedExpr<T> {
    /// Create an owned expression from an arena and a root pointer.
    ///
    /// # Safety
    ///
    /// `root` must point to a value allocated in `arena` and must be valid
    /// for the lifetime of `arena`.
    pub unsafe fn new(arena: Box<Arena>, root: *mut T) -> Self {
        Self {
            arena,
            root,
            _marker: PhantomData,
        }
    }

    /// Access the root expression.
    pub fn root(&self) -> &T {
        // SAFETY: root is valid as long as arena is alive.
        unsafe { &*self.root }
    }
}

unsafe impl<T: Send> Send for OwnedExpr<T> {}
unsafe impl<T: Sync> Sync for OwnedExpr<T> {}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn allocate_with_integers() {
        let arena = Arena::new();
        let a = arena.allocate_with(|| 42);
        let b = arena.allocate_with(|| 7);
        assert_eq!(*a, 42);
        assert_eq!(*b, 7);
    }

    #[test]
    fn allocate_with_larger_than_block() {
        let arena = Arena::with_capacity(16);
        let data = [0u8; 128];
        let ptr = arena.allocate_with(|| data);
        assert_eq!(*ptr, data);
    }

    #[test]
    fn allocate_with_many_values_triggers_new_chunks() {
        let arena = Arena::with_capacity(32);
        let mut sum = 0i64;
        for i in 0..100 {
            let value = arena.allocate_with(|| i);
            sum += *value;
        }
        assert_eq!(sum, 4950);
    }

    #[test]
    fn owned_expr_keeps_arena_alive() {
        let arena = Box::new(Arena::new());
        let root = arena.allocate_with(|| 123);
        let root_ptr: *mut i32 = root;
        // SAFETY: root was allocated in arena, and arena outlives OwnedExpr.
        let owned = unsafe { OwnedExpr::new(arena, root_ptr) };
        assert_eq!(*owned.root(), 123);
    }

    #[test]
    fn copy_values_survive_arena_drop() {
        let value = {
            let arena = Arena::new();
            let ptr = arena.allocate_with(|| 42i32);
            *ptr
        };
        assert_eq!(value, 42);
    }

    proptest! {
        #[test]
        fn allocate_random_sizes(size in 1usize..10_000) {
            let arena = Arena::with_capacity(256);
            let data: Vec<u8> = (0..size).map(|i| (i % 256) as u8).collect();
            let ptr = arena.allocate_with(|| data.clone());
            prop_assert_eq!(&ptr[..], &data[..]);
        }

        #[test]
        fn allocate_many_random_values(sizes in prop::collection::vec(1usize..512, 1..50)) {
            let arena = Arena::with_capacity(256);
            let mut total = 0usize;
            for (idx, size) in sizes.iter().enumerate() {
                let expected: Vec<u8> = (0..*size).map(|i| (i.wrapping_add(idx)) as u8).collect();
                let ptr = arena.allocate_with(|| expected.clone());
                prop_assert_eq!(&ptr[..], &expected[..]);
                total += size;
            }
            prop_assert!(total > 0);
        }
    }
}
