// oCAS C++ RAII wrapper.
//
// This header provides a thin RAII layer over the C API in ocas.h. Include
// it after ocas.h:
//
//     #include "ocas.h"
//     #include "ocas.hpp"
//
// All resources owned by ocas::Expression are released automatically when
// the object is destroyed or reassigned. Errors from the C API are
// translated into ocas::Error exceptions.

#ifndef OCAS_HPP
#define OCAS_HPP

#include <ocas.h>
#include <stdexcept>
#include <string>
#include <utility>

namespace ocas {

/// Exception thrown when an oCAS C API call fails.
class Error : public std::runtime_error {
public:
    Error(const std::string& what) : std::runtime_error(what) {}
};

/// RAII wrapper around an opaque expression handle.
class Expression {
public:
    /// Parse a string into an expression. Throws [`Error`] on parse failure.
    explicit Expression(const std::string& input) {
        int err = 0;
        handle_ = ::ocas_expr_parse(input.c_str(), &err);
        if (handle_ == nullptr) {
            throw Error(error_message());
        }
    }

    /// Take ownership of an existing handle (may be null).
    explicit Expression(::ocas_OcasExpr* handle) noexcept : handle_(handle) {}

    /// Copy constructor: clones the underlying expression.
    Expression(const Expression& other) {
        int err = 0;
        handle_ = ::ocas_expr_clone(other.handle_, &err);
        if (handle_ == nullptr) {
            throw Error(error_message());
        }
    }

    /// Copy assignment.
    Expression& operator=(const Expression& other) {
        if (this != &other) {
            Expression tmp(other);
            swap(tmp);
        }
        return *this;
    }

    /// Move constructor.
    Expression(Expression&& other) noexcept : handle_(other.handle_) {
        other.handle_ = nullptr;
    }

    /// Move assignment.
    Expression& operator=(Expression&& other) noexcept {
        if (this != &other) {
            free();
            handle_ = other.handle_;
            other.handle_ = nullptr;
        }
        return *this;
    }

    /// Destructor: releases the underlying handle.
    ~Expression() { free(); }

    /// Render the expression to a string.
    std::string to_string() const {
        int err = 0;
        char* s = ::ocas_expr_to_string(handle_, &err);
        if (s == nullptr) {
            throw Error(error_message());
        }
        std::string result(s);
        ::ocas_string_free(s);
        return result;
    }

    /// Differentiate with respect to `var`.
    Expression diff(const std::string& var) const {
        int err = 0;
        ::ocas_OcasExpr* result = ::ocas_expr_diff(handle_, var.c_str(), &err);
        if (result == nullptr) {
            throw Error(error_message());
        }
        return Expression(result);
    }

    /// Integrate with respect to `var`.
    Expression integrate(const std::string& var) const {
        int err = 0;
        ::ocas_OcasExpr* result = ::ocas_expr_integrate(handle_, var.c_str(), &err);
        if (result == nullptr) {
            throw Error(error_message());
        }
        return Expression(result);
    }

    /// Simplify using the default rule set.
    Expression simplify() const {
        int err = 0;
        ::ocas_OcasExpr* result = ::ocas_expr_simplify(handle_, &err);
        if (result == nullptr) {
            throw Error(error_message());
        }
        return Expression(result);
    }

    /// Substitute every occurrence of `var` with `replacement`.
    Expression substitute(const std::string& var, const Expression& replacement) const {
        int err = 0;
        ::ocas_OcasExpr* result =
            ::ocas_expr_substitute(handle_, var.c_str(), replacement.raw(), &err);
        if (result == nullptr) {
            throw Error(error_message());
        }
        return Expression(result);
    }

    /// Access the raw opaque handle (non-owning).
    ::ocas_OcasExpr* raw() const noexcept { return handle_; }

    /// Swap two expressions.
    void swap(Expression& other) noexcept { std::swap(handle_, other.handle_); }

private:
    ::ocas_OcasExpr* handle_;

    void free() noexcept {
        if (handle_ != nullptr) {
            ::ocas_expr_free(handle_);
            handle_ = nullptr;
        }
    }

    static std::string error_message() {
        const char* msg = ::ocas_error_last_message();
        return msg != nullptr ? std::string(msg) : std::string("unknown oCAS error");
    }
};

/// Number-theory helpers over the `ocas_ntheory_*` C API. Integers of
/// arbitrary size are passed as decimal strings; string results are
/// released automatically.
namespace ntheory {
namespace detail {
    inline std::string last_error() {
        const char* msg = ::ocas_error_last_message();
        return msg != nullptr ? std::string(msg) : std::string("unknown oCAS error");
    }

    inline std::string take_string(char* s) {
        if (s == nullptr) {
            throw Error(last_error());
        }
        std::string out(s);
        ::ocas_string_free(s);
        return out;
    }
}  // namespace detail

/// Factor `|n|` into primes as `"p1:e1,p2:e2,..."` (`"-1:1"` first when
/// negative).
inline std::string factorint(const std::string& n) {
    int err = 0;
    return detail::take_string(::ocas_ntheory_factorint(n.c_str(), &err));
}

/// BPSW probable-prime test.
inline bool isprime(const std::string& n) {
    int err = 0;
    int r = ::ocas_ntheory_isprime(n.c_str(), &err);
    if (r < 0) {
        throw Error(detail::last_error());
    }
    return r != 0;
}

/// Smallest prime strictly greater than `n`.
inline std::string nextprime(const std::string& n) {
    int err = 0;
    return detail::take_string(::ocas_ntheory_nextprime(n.c_str(), &err));
}

/// Solve `base^x ≡ target (mod p)`; throws when no logarithm exists.
inline std::string discrete_log(const std::string& p,
                                const std::string& base,
                                const std::string& target) {
    int err = 0;
    return detail::take_string(
        ::ocas_ntheory_discrete_log(p.c_str(), base.c_str(), target.c_str(), &err));
}

/// Chinese remainder theorem over comma-separated lists; returns `"r,m"`.
inline std::string crt(const std::string& moduli, const std::string& residues) {
    int err = 0;
    return detail::take_string(::ocas_ntheory_crt(moduli.c_str(), residues.c_str(), &err));
}

/// The Jacobi symbol `(a / n)` for odd positive `n`.
inline int jacobi(const std::string& a, const std::string& n) {
    int err = 0;
    int r = ::ocas_ntheory_jacobi(a.c_str(), n.c_str(), &err);
    if (r == -2) {
        throw Error(detail::last_error());
    }
    return r;
}

/// Euler's totient `φ(n)`.
inline std::string totient(const std::string& n) {
    int err = 0;
    return detail::take_string(::ocas_ntheory_totient(n.c_str(), &err));
}

/// The Möbius function `μ(n)`.
inline int mobius(const std::string& n) {
    int err = 0;
    int r = ::ocas_ntheory_mobius(n.c_str(), &err);
    if (r == -2) {
        throw Error(detail::last_error());
    }
    return r;
}

/// Number of positive divisors `τ(n)`.
inline std::string divisor_count(const std::string& n) {
    int err = 0;
    return detail::take_string(::ocas_ntheory_divisor_count(n.c_str(), &err));
}

/// Sum of `k`-th powers of the positive divisors `σ_k(n)`.
inline std::string divisor_sigma(const std::string& n, uint32_t k = 1) {
    int err = 0;
    return detail::take_string(::ocas_ntheory_divisor_sigma(n.c_str(), k, &err));
}

/// Liouville's function `λ(n)`.
inline int liouville(const std::string& n) {
    int err = 0;
    int r = ::ocas_ntheory_liouville(n.c_str(), &err);
    if (r == -2) {
        throw Error(detail::last_error());
    }
    return r;
}
}  // namespace ntheory

/// Tensor algebra helpers (0.22.0) over the `ocas_tensor_*` C API.
namespace tensor {

/// Canonicalise a tensor expression via graph isomorphism.
/// `specs` is a comma-separated `name:sym` string, e.g. `"T:none,U:antisymmetric"`.
/// `groups` (optional) is a comma-separated `label:group` string.
inline std::string canonicalize(const std::string& expr,
                                const std::string& specs,
                                const std::string& groups = "") {
    int err = 0;
    char* s = ::ocas_tensor_canonicalize(
        expr.c_str(),
        specs.c_str(),
        groups.empty() ? nullptr : groups.c_str(),
        &err);
    return ntheory::detail::take_string(s);
}

/// Apply a Young tableau projector. `tableau` is a comma-separated list
/// of row lengths, e.g. `"2,1"` for □□/□.
inline std::string young_project(const std::string& expr,
                                 const std::string& tableau) {
    int err = 0;
    char* s = ::ocas_young_project(expr.c_str(), tableau.c_str(), &err);
    return ntheory::detail::take_string(s);
}

/// Refresh dummy indices in a tensor expression.
inline std::string refresh_dummies(const std::string& expr,
                                   const std::string& specs) {
    int err = 0;
    char* s = ::ocas_tensor_refresh_dummies(expr.c_str(), specs.c_str(), &err);
    return ntheory::detail::take_string(s);
}

}  // namespace tensor

}  // namespace ocas

#endif  // OCAS_HPP
