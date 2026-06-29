#include <stdio.h>
#include <stdlib.h>
#include "../include/ocas.h"

int main(void) {
    const char *version = ocas_version();
    printf("oCAS version: %s\n", version);

    struct ocas_OcasArena *arena = ocas_arena_new();
    if (arena == NULL) {
        fprintf(stderr, "failed to create arena: %s\n", ocas_error_last_message());
        return 1;
    }

    int64_t value = 0;
    int rc = ocas_arena_alloc_i64(arena, 42, &value);
    if (rc != ocas_OCAS_OK) {
        fprintf(stderr, "failed to allocate in arena: %s\n", ocas_error_last_message());
        ocas_arena_free(arena);
        return 1;
    }

    printf("allocated value: %lld\n", (long long)value);

    struct ocas_OcasArena *null_arena = NULL;
    int null_rc = ocas_arena_alloc_i64(null_arena, 0, &value);
    if (null_rc == ocas_OCAS_OK) {
        fprintf(stderr, "expected null arena to return an error\n");
        ocas_arena_free(arena);
        return 1;
    }
    printf("null arena error: %s\n", ocas_error_last_message());
    ocas_error_clear();

    ocas_arena_free(arena);
    printf("basic C example completed successfully\n");
    return 0;
}
