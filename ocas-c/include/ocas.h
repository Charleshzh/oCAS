/* oCAS C API */

/* Generated with cbindgen:0.28.0 */

/* Warning: this file is auto-generated. Do not modify. */

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace ocas {
#endif  // __cplusplus

/**
 * Success error code.
 */
#define ocas_OCAS_OK 0

/**
 * A null pointer was passed where a non-null pointer was required.
 */
#define ocas_OCAS_ERROR_NULL_POINTER 1

/**
 * An operation failed inside the oCAS runtime.
 */
#define ocas_OCAS_ERROR_RUNTIME 2

/**
 * Opaque arena handle.
 */
typedef struct ocas_OcasArena {
  uint8_t _private[0];
} ocas_OcasArena;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Return the oCAS version string.
 *
 * The returned pointer is valid for the lifetime of the program and must not
 * be freed by the caller.
 */
const char *ocas_version(void);

/**
 * Return the message for the last error on the calling thread, or `NULL` if
 * no error has occurred.
 *
 * The returned string is owned by the library and must not be freed or
 * modified by the caller. It remains valid until the next call that sets an
 * error on the same thread or until `ocas_error_clear` is called.
 */
const char *ocas_error_last_message(void);

/**
 * Clear the last error on the calling thread.
 */
void ocas_error_clear(void);

/**
 * Create a new arena and return an opaque pointer to it.
 *
 * Returns `NULL` if allocation fails. Use `ocas_error_last_message` to
 * retrieve the error message.
 */
struct ocas_OcasArena *ocas_arena_new(void);

/**
 * Free an arena previously created with `ocas_arena_new`.
 *
 * Passing `NULL` is a no-op.
 */
void ocas_arena_free(struct ocas_OcasArena *arena);

/**
 * Allocate a single `i64` in the arena and return its value.
 *
 * This is a trivial demonstration of the arena lifetime model. Returns
 * `OCAS_ERROR_NULL_POINTER` if `arena` is null; otherwise returns `OCAS_OK`.
 *
 * # Safety
 *
 * `arena` must be a non-null pointer returned by `ocas_arena_new`. `out`
 * must be a valid, non-null pointer to writable memory.
 */
int ocas_arena_alloc_i64(struct ocas_OcasArena *arena, int64_t value, int64_t *out);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

#ifdef __cplusplus
}  // namespace ocas
#endif  // __cplusplus
