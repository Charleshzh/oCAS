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
///
/// # Example
///
/// ```
/// use ocas_core::arena::Arena;
///
/// let arena = Arena::new();
/// let value = arena.allocate_with(|| 42);
/// assert_eq!(*value, 42);
/// ```
pub struct Arena {
    chunks: RefCell<Vec<Chunk>>,
    block_size: usize,
}

impl Arena {
    /// Create a new arena with the default block size.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let n = arena.allocate_with(|| 7);
    /// assert_eq!(*n, 7);
    /// ```
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
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let value = arena.allocate_with(|| "hello");
    /// assert_eq!(*value, "hello");
    /// ```
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

    /// Allocate a contiguous slice of `T` values in the arena.
    ///
    /// The returned slice is tied to the arena lifetime. Because the arena does
    /// not run destructors, `T` must be `Copy` so that dropping the arena does
    /// not leak resources owned by the slice elements.
    ///
    /// # Panics
    ///
    /// Panics if `T` has zero size or if the total allocation size overflows.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let slice = arena.allocate_slice(&[1, 2, 3, 4, 5]);
    /// assert_eq!(slice, &[1, 2, 3, 4, 5]);
    /// ```
    pub fn allocate_slice<T: Copy>(&self, values: &[T]) -> &[T] {
        if values.is_empty() {
            return &[];
        }

        let layout = Layout::from_size_align(mem::size_of_val(values), mem::align_of::<T>())
            .expect("invalid slice layout");
        assert!(
            layout.size() > 0,
            "cannot allocate zero-sized types in Arena"
        );

        let ptr = self.alloc_raw(layout);

        // SAFETY: `ptr` is non-null, properly aligned, and points to a block
        // large enough to hold `values.len()` elements of type `T`.
        unsafe {
            let typed = ptr.as_ptr().cast::<T>();
            std::ptr::copy_nonoverlapping(values.as_ptr(), typed, values.len());
            std::slice::from_raw_parts(typed, values.len())
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

        // Need a new chunk. Use at least the requested size and alignment so the
        // first allocation in the chunk is correctly aligned.
        let size = layout.size().max(self.block_size);
        let align = layout.align();
        let mut new_chunk = Chunk::new(size, align);
        let ptr = new_chunk
            .try_alloc(layout)
            .expect("new chunk should fit any layout up to its size");
        chunks.push(new_chunk);
        ptr
    }

    /// Reset the arena, invalidating all previously allocated values.
    ///
    /// The first chunk is kept and reused; additional chunks are released.
    /// This makes repeated build–reset cycles allocation-free in the steady
    /// state, which is the basis of the workspace pool in `ocas-atom`.
    ///
    /// # Safety contract (enforced by convention)
    ///
    /// Any reference returned by [`allocate_with`](Arena::allocate_with) or
    /// [`allocate_slice`](Arena::allocate_slice) before the reset **must not**
    /// be used afterwards — the memory may be handed out again for different
    /// values. Callers must treat reset as the end of a generation.
    ///
    /// # Example
    ///
    /// ```
    /// use ocas_core::arena::Arena;
    ///
    /// let arena = Arena::new();
    /// let _ = arena.allocate_with(|| 1);
    /// arena.reset();
    /// let value = arena.allocate_with(|| 2);
    /// assert_eq!(*value, 2);
    /// ```
    pub fn reset(&self) {
        let mut chunks = self.chunks.borrow_mut();
        // Keep the first chunk (reusable block), release the rest.
        chunks.truncate(1);
        if let Some(first) = chunks.first_mut() {
            first.offset = 0;
        }
    }

    /// Return the number of chunks currently held by the arena.
    pub fn chunk_count(&self) -> usize {
        self.chunks.borrow().len()
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
    align: usize,
    offset: usize,
}

impl Chunk {
    fn new(size: usize, align: usize) -> Self {
        let layout = Layout::from_size_align(size, align).expect("invalid chunk layout");
        // SAFETY: layout is non-zero and properly aligned.
        let memory = unsafe { NonNull::new_unchecked(alloc::alloc(layout)) };
        Self {
            memory,
            size,
            align,
            offset: 0,
        }
    }

    fn try_alloc(&mut self, layout: Layout) -> Option<NonNull<u8>> {
        // The chunk's base pointer is only guaranteed to be aligned to
        // `self.align`; refuse requests needing stricter alignment so the
        // caller falls through to a fresh, suitably-aligned chunk.
        if layout.align() > self.align {
            return None;
        }
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
        // Deallocation requires the same size and an alignment that is at least
        // as large as the original allocation. The chunk was allocated with the
        // maximum alignment requested by any layout served from this chunk, so
        // that alignment is stored alongside the chunk.
        let layout = Layout::from_size_align(self.size, self.align).expect("invalid chunk layout");
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

    mod reset {
        use super::*;

        #[test]
        fn reset_keeps_first_chunk() {
            let arena = Arena::with_capacity(64);
            // Force several chunks (64-byte blocks, 8-byte values).
            for i in 0..100u64 {
                arena.allocate_with(|| i);
            }
            assert!(arena.chunk_count() > 1);
            arena.reset();
            assert_eq!(arena.chunk_count(), 1);
        }

        #[test]
        fn reset_allows_reuse() {
            let arena = Arena::new();
            let _ = arena.allocate_with(|| 1u64);
            arena.reset();
            let value = arena.allocate_with(|| 2u64);
            assert_eq!(*value, 2);
        }

        #[test]
        fn reset_reuse_steady_state_allocates_nothing() {
            let arena = Arena::with_capacity(4096);
            for round in 0..1000 {
                for i in 0..50u64 {
                    let v = arena.allocate_with(|| i + round);
                    assert_eq!(*v, i + round);
                }
                arena.reset();
            }
            // Steady state: still a single chunk after 1000 generations.
            assert_eq!(arena.chunk_count(), 1);
        }

        #[test]
        fn overaligned_allocation_gets_own_chunk() {
            let arena = Arena::new();
            let _ = arena.allocate_with(|| 1u8);
            #[repr(align(64))]
            #[derive(Copy, Clone)]
            struct Wide(u64);
            let w = arena.allocate_with(|| Wide(7));
            assert_eq!(w as *const Wide as usize % 64, 0);
            assert_eq!(w.0, 7);
        }
    }

    mod simple {
        use super::*;

        #[test]
        fn allocate_single_integer() {
            let arena = Arena::new();
            let value = arena.allocate_with(|| 42);
            assert_eq!(*value, 42);
        }

        #[test]
        fn allocate_two_integers() {
            let arena = Arena::new();
            let a = arena.allocate_with(|| 1);
            let b = arena.allocate_with(|| 2);
            assert_eq!(*a, 1);
            assert_eq!(*b, 2);
        }

        #[test]
        fn allocate_empty_slice() {
            let arena = Arena::new();
            let slice: &[i32] = arena.allocate_slice(&[]);
            assert!(slice.is_empty());
        }

        #[test]
        fn allocate_small_slice() {
            let arena = Arena::new();
            let data = [10, 20, 30];
            let slice = arena.allocate_slice(&data);
            assert_eq!(slice, &data[..]);
        }

        #[test]
        fn arena_default_matches_new() {
            let arena: Arena = Default::default();
            let value = arena.allocate_with(|| "x");
            assert_eq!(*value, "x");
        }
    }

    mod medium {
        use super::*;

        #[test]
        fn allocate_larger_than_block() {
            let arena = Arena::with_capacity(16);
            let data = [0u8; 128];
            let ptr = arena.allocate_with(|| data);
            assert_eq!(*ptr, data);
        }

        #[test]
        fn allocate_slice_larger_than_block() {
            let arena = Arena::with_capacity(16);
            let values: Vec<u8> = (0..=255).collect();
            let slice = arena.allocate_slice(&values);
            assert_eq!(slice, &values[..]);
        }

        #[test]
        fn multiple_chunks_for_many_values() {
            let arena = Arena::with_capacity(32);
            let mut sum = 0i64;
            for i in 0..100 {
                let value = arena.allocate_with(|| i);
                sum += *value;
            }
            assert_eq!(sum, 4950);
        }

        #[test]
        fn multiple_chunks_for_many_slices() {
            let arena = Arena::with_capacity(64);
            let mut total = 0i64;
            for i in 0..50 {
                let values: Vec<i64> = (0..10).map(|j| i * 10 + j).collect();
                let slice = arena.allocate_slice(&values);
                total += slice.iter().sum::<i64>();
            }
            assert_eq!(total, 124_750);
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
    }

    mod complex {
        use super::*;

        #[test]
        fn copy_values_survive_arena_drop() {
            let value = {
                let arena = Arena::new();
                let ptr = arena.allocate_with(|| 42i32);
                *ptr
            };
            assert_eq!(value, 42);
        }

        #[test]
        fn alignment_of_large_type() {
            #[derive(Clone, Copy)]
            #[repr(C, align(64))]
            struct BigAlign(u64);

            // The first allocation in a fresh chunk must respect the requested
            // alignment. Use a single-element slice so the layout alignment is
            // dominated by BigAlign.
            let arena = Arena::with_capacity(4096);
            let values = [BigAlign(7)];
            let slice = arena.allocate_slice(&values);
            assert!((slice.as_ptr() as usize).is_multiple_of(64));
            assert_eq!(slice[0].0, 7);
        }

        #[test]
        fn alignment_of_single_value() {
            #[derive(Clone, Copy)]
            #[repr(C, align(64))]
            struct BigAlign(u64);

            let arena = Arena::with_capacity(4096);
            let value = arena.allocate_with(|| BigAlign(7));
            assert_eq!((value as *const BigAlign) as usize % 64, 0);
            assert_eq!(value.0, 7);
        }

        #[test]
        #[should_panic(expected = "cannot allocate zero-sized types in Arena")]
        fn zero_sized_type_panics() {
            let arena = Arena::new();
            let _: &mut () = arena.allocate_with(|| ());
        }

        #[test]
        #[should_panic(expected = "cannot allocate zero-sized types in Arena")]
        fn zero_sized_slice_panics() {
            let arena = Arena::new();
            let _: &[()] = arena.allocate_slice(&[()]);
        }

        #[test]
        fn owned_expr_is_send_sync() {
            fn assert_send_sync<T: Send + Sync>() {}
            assert_send_sync::<OwnedExpr<u8>>();
        }
    }

    mod extreme {
        use super::*;

        #[test]
        fn stress_mixed_allocations() {
            let arena = Arena::with_capacity(256);
            let mut total = 0usize;
            for size in (1usize..=1000).step_by(7) {
                let data: Vec<u8> = (0..size).map(|i| (i % 256) as u8).collect();
                let ptr = arena.allocate_with(|| data.clone());
                total += ptr.iter().map(|&x| x as usize).sum::<usize>();
            }
            assert!(total > 0);
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

            #[test]
            fn slice_roundtrip(values in prop::collection::vec(0i32..100, 0..512)) {
                let arena = Arena::with_capacity(256);
                let slice = arena.allocate_slice(&values);
                prop_assert_eq!(slice, &values[..]);
            }
        }
    }
}
