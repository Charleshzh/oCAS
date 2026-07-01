// Example: oCAS C++ RAII wrapper.
//
// Compile and link against libocas_c plus this file. Demonstrates
// automatic resource management, differentiation, integration, and
// substitution through the ocas::Expression class.

#include <iostream>
#include <ocas.hpp>

int main() {
    try {
        ocas::Expression expr("x^3");
        std::cout << "expr: " << expr.to_string() << '\n';

        // d/dx(x^3)
        ocas::Expression deriv = expr.diff("x");
        std::cout << "d/dx(x^3) = " << deriv.to_string() << '\n';

        // Second derivative.
        ocas::Expression second = deriv.diff("x");
        std::cout << "d^2/dx^2(x^3) = " << second.to_string() << '\n';

        // Substitute x -> 2 in the second derivative: 6*2 = 12.
        ocas::Expression replacement("2");
        ocas::Expression evaluated = second.substitute("x", replacement);
        std::cout << "d^2/dx^2(x^3) at x=2 = " << evaluated.to_string() << '\n';

        // Copy and move semantics.
        ocas::Expression copy = expr;  // copy
        ocas::Expression moved = std::move(copy);  // move
        std::cout << "moved: " << moved.to_string() << '\n';

    } catch (const ocas::Error& e) {
        std::cerr << "oCAS error: " << e.what() << '\n';
        return 1;
    }

    std::cout << "C++ example completed successfully\n";
    return 0;
}
