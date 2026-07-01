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

}  // namespace ocas

#endif  // OCAS_HPP
