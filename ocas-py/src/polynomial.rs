//! Python `Polynomial` class — dense univariate polynomials over ℤ, ℚ, or GF(p).
//!
//! Wraps [`ocas_poly::DenseUnivariatePolynomial`] with an enum-erasure
//! strategy so that a single Python class supports three coefficient domains.

use ocas_domain::{FiniteField, Integer, IntegerDomain, Rational, RationalDomain};
use ocas_poly::DenseUnivariatePolynomial;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

use crate::domain::DomainKind;

/// Type-erased polynomial over one of the three supported domains.
#[derive(Clone)]
pub(crate) enum PolyErased {
    Int(DenseUnivariatePolynomial<IntegerDomain>),
    Rat(DenseUnivariatePolynomial<RationalDomain>),
    Fq(DenseUnivariatePolynomial<FiniteField>),
}

/// A dense univariate polynomial.
///
/// The coefficient domain is selected by the `domain` argument: one of
/// the strings `"integer"` (default), `"rational"`, or a `FiniteField`
/// instance.
///
/// ```python
/// from ocas import Polynomial
///
/// # x^2 + 2x + 1 over the integers
/// p = Polynomial([1, 2, 1])
/// print(p.degree())       # 2
/// print(p.eval(2))        # 9
/// print((p * p).coeffs()) # [1, 4, 6, 4, 1]
///
/// # over GF(5)
/// from ocas import FiniteField
/// q = Polynomial([1, 2, 1], domain=FiniteField(5))
/// ```
#[pyclass(name = "Polynomial", skip_from_py_object)]
#[derive(Clone)]
pub struct PyPolynomial {
    pub(crate) inner: PolyErased,
}

/// A single (polynomial, multiplicity) factor returned by factorization.
#[pyclass(name = "PolynomialFactor", skip_from_py_object)]
pub struct PyPolynomialFactor {
    #[pyo3(get)]
    pub factor: PyPolynomial,
    #[pyo3(get)]
    pub multiplicity: usize,
}

/// Extract integer-domain coefficients from a Python iterable of ints.
fn extract_int_coeffs(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Integer>> {
    let ints: Vec<i64> = obj
        .extract()
        .map_err(|_| PyTypeError::new_err("integer coefficients must be ints"))?;
    Ok(ints.into_iter().map(Integer::from).collect())
}

/// Extract rational-domain coefficients. Each element is either an int
/// (denominator 1) or a `(numerator, denominator)` tuple; the whole list
/// must be uniformly one form.
fn extract_rat_coeffs(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Rational>> {
    if let Ok(ints) = obj.extract::<Vec<i64>>() {
        Ok(ints.into_iter().map(|n| Rational::new(n, 1)).collect())
    } else if let Ok(pairs) = obj.extract::<Vec<(i64, i64)>>() {
        pairs
            .into_iter()
            .map(|(num, den)| {
                if den == 0 {
                    Err(PyValueError::new_err("rational denominator cannot be zero"))
                } else {
                    Ok(Rational::new(num, den))
                }
            })
            .collect()
    } else {
        Err(PyTypeError::new_err(
            "rational coefficients must be ints or (num, denom) tuples",
        ))
    }
}

/// Build a `PyPolynomial` from a Python iterable of coefficients and a domain kind.
pub(crate) fn build_polynomial(
    coeffs: &Bound<'_, PyAny>,
    domain: &DomainKind,
) -> PyResult<PyPolynomial> {
    let inner = match domain {
        DomainKind::Integer => {
            let c = extract_int_coeffs(coeffs)?;
            PolyErased::Int(DenseUnivariatePolynomial::from_coeffs(IntegerDomain, c))
        }
        DomainKind::Rational => {
            let c = extract_rat_coeffs(coeffs)?;
            PolyErased::Rat(DenseUnivariatePolynomial::from_coeffs(RationalDomain, c))
        }
        DomainKind::FiniteField(p) => {
            let field = FiniteField::new(p.clone());
            let ints: Vec<i64> = coeffs
                .extract()
                .map_err(|_| PyTypeError::new_err("finite-field coefficients must be ints"))?;
            let c: Vec<_> = ints.into_iter().map(|v| field.element(v)).collect();
            PolyErased::Fq(DenseUnivariatePolynomial::from_coeffs(field, c))
        }
    };
    Ok(PyPolynomial { inner })
}

#[pymethods]
impl PyPolynomial {
    /// Create a polynomial from coefficients (constant term first).
    ///
    /// `domain` selects the coefficient ring: `"integer"` (default),
    /// `"rational"`, or a `FiniteField` instance.
    #[new]
    #[pyo3(signature = (coeffs, domain=None))]
    fn new(coeffs: &Bound<'_, PyAny>, domain: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let kind = match domain {
            Some(d) => DomainKind::from_py(d)?,
            None => DomainKind::Integer,
        };
        build_polynomial(coeffs, &kind)
    }

    /// Return the coefficients as a list of decimal strings (constant term
    /// first). String form preserves arbitrary precision across the
    /// gmp/non-gmp builds; wrap each entry in `int(...)` to obtain a Python
    /// integer. Rational entries are rendered as `n/d`.
    fn coeffs(&self) -> Vec<String> {
        match &self.inner {
            PolyErased::Int(p) => p.coeffs().iter().map(|c| c.to_string()).collect(),
            PolyErased::Rat(p) => p.coeffs().iter().map(|c| c.to_string()).collect(),
            PolyErased::Fq(p) => p.coeffs().iter().map(|c| c.value().to_string()).collect(),
        }
    }

    /// Return the degree, or `None` for the zero polynomial.
    fn degree(&self) -> Option<usize> {
        match &self.inner {
            PolyErased::Int(p) => p.degree(),
            PolyErased::Rat(p) => p.degree(),
            PolyErased::Fq(p) => p.degree(),
        }
    }

    /// Return the number of stored coefficients.
    fn len(&self) -> usize {
        match &self.inner {
            PolyErased::Int(p) => p.coeffs().len(),
            PolyErased::Rat(p) => p.coeffs().len(),
            PolyErased::Fq(p) => p.coeffs().len(),
        }
    }

    /// Return `True` if this is the zero polynomial.
    fn is_zero(&self) -> bool {
        self.len() == 0
    }

    /// Evaluate the polynomial at `x` and return the result as a decimal
    /// string (rational results are rendered as `n/d`).
    ///
    /// For integer/finite-field domains, `x` is an int. For the rational
    /// domain, `x` may be an int or a `(num, denom)` tuple.
    fn eval(&self, x: &Bound<'_, PyAny>) -> PyResult<String> {
        match &self.inner {
            PolyErased::Int(p) => {
                let v = x
                    .extract::<i64>()
                    .map_err(|_| PyTypeError::new_err("x must be an int"))?;
                Ok(p.eval(&Integer::from(v)).to_string())
            }
            PolyErased::Rat(p) => {
                let v = if let Ok(n) = x.extract::<i64>() {
                    Rational::new(n, 1)
                } else if let Ok((num, den)) = x.extract::<(i64, i64)>() {
                    Rational::new(num, den)
                } else {
                    return Err(PyTypeError::new_err(
                        "x must be an int or (num, denom) tuple",
                    ));
                };
                Ok(p.eval(&v).to_string())
            }
            PolyErased::Fq(p) => {
                let field = p.domain();
                let v = x
                    .extract::<i64>()
                    .map_err(|_| PyTypeError::new_err("x must be an int"))?;
                Ok(p.eval(&field.element(v)).value().to_string())
            }
        }
    }

    /// Return the formal derivative.
    fn derivative(&self) -> PyPolynomial {
        match &self.inner {
            PolyErased::Int(p) => PyPolynomial {
                inner: PolyErased::Int(p.derivative()),
            },
            PolyErased::Rat(p) => PyPolynomial {
                inner: PolyErased::Rat(p.derivative()),
            },
            PolyErased::Fq(p) => PyPolynomial {
                inner: PolyErased::Fq(p.derivative()),
            },
        }
    }

    /// Return the formal integral with constant term zero.
    fn integral(&self) -> PyPolynomial {
        match &self.inner {
            PolyErased::Int(p) => PyPolynomial {
                inner: PolyErased::Int(p.integral()),
            },
            PolyErased::Rat(p) => PyPolynomial {
                inner: PolyErased::Rat(p.integral()),
            },
            PolyErased::Fq(p) => PyPolynomial {
                inner: PolyErased::Fq(p.integral()),
            },
        }
    }

    /// Return the primitive part (content stripped) for integer polynomials.
    fn primitive_part(&self) -> PyResult<PyPolynomial> {
        match &self.inner {
            PolyErased::Int(p) => Ok(PyPolynomial {
                inner: PolyErased::Int(p.primitive_part()),
            }),
            _ => Err(PyValueError::new_err(
                "primitive_part is only defined over the integers",
            )),
        }
    }

    /// Return the square-free factorization as a list of `(factor, multiplicity)`.
    fn square_free_factorization(&self) -> PyResult<Vec<PyPolynomialFactor>> {
        let factors: Vec<_> = match &self.inner {
            PolyErased::Int(p) => p
                .square_free_factorization()
                .into_iter()
                .map(|(f, m)| PyPolynomialFactor {
                    factor: PyPolynomial {
                        inner: PolyErased::Int(f),
                    },
                    multiplicity: m,
                })
                .collect(),
            PolyErased::Rat(p) => p
                .square_free_factorization()
                .into_iter()
                .map(|(f, m)| PyPolynomialFactor {
                    factor: PyPolynomial {
                        inner: PolyErased::Rat(f),
                    },
                    multiplicity: m,
                })
                .collect(),
            PolyErased::Fq(p) => p
                .square_free_factorization()
                .into_iter()
                .map(|(f, m)| PyPolynomialFactor {
                    factor: PyPolynomial {
                        inner: PolyErased::Fq(f),
                    },
                    multiplicity: m,
                })
                .collect(),
        };
        Ok(factors)
    }

    /// Return `True` if the polynomial has no repeated factors.
    fn is_square_free(&self) -> bool {
        match &self.inner {
            PolyErased::Int(p) => p.is_square_free(),
            PolyErased::Rat(p) => p.is_square_free(),
            PolyErased::Fq(p) => p.is_square_free(),
        }
    }

    /// Return the greatest common divisor with `other`.
    ///
    /// Both polynomials must share the same coefficient domain.
    fn gcd(&self, other: &PyPolynomial) -> PyResult<PyPolynomial> {
        match (&self.inner, &other.inner) {
            (PolyErased::Int(a), PolyErased::Int(b)) => Ok(PyPolynomial {
                inner: PolyErased::Int(a.gcd(b)),
            }),
            (PolyErased::Rat(a), PolyErased::Rat(b)) => Ok(PyPolynomial {
                inner: PolyErased::Rat(a.gcd(b)),
            }),
            (PolyErased::Fq(a), PolyErased::Fq(b)) => Ok(PyPolynomial {
                inner: PolyErased::Fq(a.gcd(b)),
            }),
            _ => Err(PyTypeError::new_err(
                "gcd requires both polynomials to share the same coefficient domain",
            )),
        }
    }

    /// Divide by `other`, returning `(quotient, remainder)`, or `None` if
    /// `other` is zero.
    fn div_rem(&self, other: &PyPolynomial) -> PyResult<Option<(PyPolynomial, PyPolynomial)>> {
        match (&self.inner, &other.inner) {
            (PolyErased::Int(a), PolyErased::Int(b)) => Ok(a.div_rem(b).map(|(q, r)| {
                (
                    PyPolynomial {
                        inner: PolyErased::Int(q),
                    },
                    PyPolynomial {
                        inner: PolyErased::Int(r),
                    },
                )
            })),
            (PolyErased::Rat(a), PolyErased::Rat(b)) => Ok(a.div_rem(b).map(|(q, r)| {
                (
                    PyPolynomial {
                        inner: PolyErased::Rat(q),
                    },
                    PyPolynomial {
                        inner: PolyErased::Rat(r),
                    },
                )
            })),
            (PolyErased::Fq(a), PolyErased::Fq(b)) => Ok(a.div_rem(b).map(|(q, r)| {
                (
                    PyPolynomial {
                        inner: PolyErased::Fq(q),
                    },
                    PyPolynomial {
                        inner: PolyErased::Fq(r),
                    },
                )
            })),
            _ => Err(PyTypeError::new_err(
                "div_rem requires both polynomials to share the same coefficient domain",
            )),
        }
    }

    /// Return `self + other`.
    fn __add__(&self, other: &PyPolynomial) -> PyResult<PyPolynomial> {
        match (&self.inner, &other.inner) {
            (PolyErased::Int(a), PolyErased::Int(b)) => Ok(PyPolynomial {
                inner: PolyErased::Int(a.add(b)),
            }),
            (PolyErased::Rat(a), PolyErased::Rat(b)) => Ok(PyPolynomial {
                inner: PolyErased::Rat(a.add(b)),
            }),
            (PolyErased::Fq(a), PolyErased::Fq(b)) => Ok(PyPolynomial {
                inner: PolyErased::Fq(a.add(b)),
            }),
            _ => Err(PyTypeError::new_err(
                "+ requires both polynomials to share the same coefficient domain",
            )),
        }
    }

    /// Return `self - other`.
    fn __sub__(&self, other: &PyPolynomial) -> PyResult<PyPolynomial> {
        match (&self.inner, &other.inner) {
            (PolyErased::Int(a), PolyErased::Int(b)) => Ok(PyPolynomial {
                inner: PolyErased::Int(a.sub(b)),
            }),
            (PolyErased::Rat(a), PolyErased::Rat(b)) => Ok(PyPolynomial {
                inner: PolyErased::Rat(a.sub(b)),
            }),
            (PolyErased::Fq(a), PolyErased::Fq(b)) => Ok(PyPolynomial {
                inner: PolyErased::Fq(a.sub(b)),
            }),
            _ => Err(PyTypeError::new_err(
                "- requires both polynomials to share the same coefficient domain",
            )),
        }
    }

    /// Return `self * other`.
    fn __mul__(&self, other: &PyPolynomial) -> PyResult<PyPolynomial> {
        match (&self.inner, &other.inner) {
            (PolyErased::Int(a), PolyErased::Int(b)) => Ok(PyPolynomial {
                inner: PolyErased::Int(a.mul(b)),
            }),
            (PolyErased::Rat(a), PolyErased::Rat(b)) => Ok(PyPolynomial {
                inner: PolyErased::Rat(a.mul(b)),
            }),
            (PolyErased::Fq(a), PolyErased::Fq(b)) => Ok(PyPolynomial {
                inner: PolyErased::Fq(a.mul(b)),
            }),
            _ => Err(PyTypeError::new_err(
                "* requires both polynomials to share the same coefficient domain",
            )),
        }
    }

    /// Return `-self`.
    fn __neg__(&self) -> PyPolynomial {
        match &self.inner {
            PolyErased::Int(p) => PyPolynomial {
                inner: PolyErased::Int(p.mul_scalar(&Integer::from(-1))),
            },
            PolyErased::Rat(p) => PyPolynomial {
                inner: PolyErased::Rat(p.mul_scalar(&Rational::new(-1, 1))),
            },
            PolyErased::Fq(p) => {
                let field = p.domain();
                PyPolynomial {
                    inner: PolyErased::Fq(p.mul_scalar(&field.element(-1))),
                }
            }
        }
    }

    /// Return `True` if the normalized coefficient vectors match.
    fn __eq__(&self, other: &PyPolynomial) -> bool {
        match (&self.inner, &other.inner) {
            (PolyErased::Int(a), PolyErased::Int(b)) => a == b,
            (PolyErased::Rat(a), PolyErased::Rat(b)) => a == b,
            (PolyErased::Fq(a), PolyErased::Fq(b)) => a == b,
            _ => false,
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            PolyErased::Int(p) => {
                format!("Polynomial([{}], 'integer')", fmt_poly_coeffs(p))
            }
            PolyErased::Rat(p) => {
                format!("Polynomial([{}], 'rational')", fmt_poly_coeffs(p))
            }
            PolyErased::Fq(p) => format!(
                "Polynomial([{}], domain=FiniteField({}))",
                fmt_poly_coeffs(p),
                p.domain().prime()
            ),
        }
    }
}

/// Format polynomial coefficients as a comma-separated list of stringified
/// values, used by `__repr__`.
fn fmt_poly_coeffs<D: ocas_domain::Domain>(p: &DenseUnivariatePolynomial<D>) -> String
where
    D::Element: std::fmt::Display,
{
    p.coeffs()
        .iter()
        .map(|c| c.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}
