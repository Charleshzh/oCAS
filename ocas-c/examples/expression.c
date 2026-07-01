// Example: expression lifecycle and calculus operations via the oCAS C API.
//
// Compile and link against libocas_c. Demonstrates parse → diff →
// to_string → free, plus substitution and error handling.

#include <stdio.h>
#include <stdlib.h>
#include "../include/ocas.h"

int main(void) {
    int err = 0;

    // Parse x^2.
    struct ocas_OcasExpr *expr = ocas_expr_parse("x^2", &err);
    if (expr == NULL) {
        fprintf(stderr, "parse failed: %s\n", ocas_error_last_message());
        return 1;
    }

    // d/dx(x^2) = 2*x.
    struct ocas_OcasExpr *deriv = ocas_expr_diff(expr, "x", &err);
    if (deriv == NULL) {
        fprintf(stderr, "diff failed: %s\n", ocas_error_last_message());
        ocas_expr_free(expr);
        return 1;
    }
    char *deriv_str = ocas_expr_to_string(deriv, &err);
    printf("d/dx(x^2) = %s\n", deriv_str);
    ocas_string_free(deriv_str);

    // Integrate: ∫ 2*x dx.
    struct ocas_OcasExpr *integ = ocas_expr_integrate(deriv, "x", &err);
    if (integ == NULL) {
        fprintf(stderr, "integrate failed: %s\n", ocas_error_last_message());
        ocas_expr_free(deriv);
        ocas_expr_free(expr);
        return 1;
    }
    char *integ_str = ocas_expr_to_string(integ, &err);
    printf("integral of 2*x dx = %s\n", integ_str);
    ocas_string_free(integ_str);

    // Substitute x -> 2 in x^2: result is 4.
    struct ocas_OcasExpr *two = ocas_expr_parse("2", &err);
    struct ocas_OcasExpr *subst = ocas_expr_substitute(expr, "x", two, &err);
    if (subst != NULL) {
        char *subst_str = ocas_expr_to_string(subst, &err);
        if (subst_str != NULL) {
            printf("x^2 with x=2 = %s\n", subst_str);
            ocas_string_free(subst_str);
        } else {
            fprintf(stderr, "to_string failed: %s\n", ocas_error_last_message());
        }
    }

    // Demonstrate NULL handle error.
    struct ocas_OcasExpr *bad = ocas_expr_diff(NULL, "x", &err);
    if (bad == NULL) {
        printf("expected error on NULL handle: %s\n", ocas_error_last_message());
        ocas_error_clear();
    }

    ocas_expr_free(subst);
    ocas_expr_free(two);
    ocas_expr_free(integ);
    ocas_expr_free(deriv);
    ocas_expr_free(expr);
    printf("expression C example completed successfully\n");
    return 0;
}
