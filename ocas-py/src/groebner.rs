//! Python bindings for Gröbner basis computation and ideal operations.

use ocas_domain::{Domain, Rational, RationalDomain};
use ocas_poly::ideal::{self, PolynomialSystemSolution};
use ocas_poly::sparse::Lex;
use ocas_poly::{
    Algorithm, GroebnerBasis, SparseMultivariatePolynomial, eliminate, groebner_basis,
};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

use crate::polynomial::PyPolynomial;

/// A multivariate polynomial over $\mathbb{Q}$, constructed from a dictionary
/// mapping exponent tuples to coefficients.
///
/// ```python
/// from ocas import MultivariatePolynomial
///
/// # x² + y² - 1 in k[x,y]
/// p = MultivariatePolynomial({(2, 0): 1, (0, 2): 1, (0, 0): -1}, n_vars=2)
/// ```
#[pyclass(name = "MultivariatePolynomial")]
pub struct PyMultivariatePolynomial {
    pub inner: SparseMultivariatePolynomial<RationalDomain, Lex>,
}

#[pymethods]
impl PyMultivariatePolynomial {
    #[new]
    fn new(terms: &Bound<'_, PyAny>, n_vars: usize) -> PyResult<Self> {
        let dict: &Bound<'_, pyo3::types::PyDict> = terms.cast().map_err(|_| {
            PyTypeError::new_err("expected a dict mapping exponent tuples to coefficients")
        })?;

        let domain = RationalDomain;
        let mut poly_terms: Vec<(Vec<usize>, Rational)> = Vec::new();

        for (key, val) in dict.iter() {
            let exp_tuple: &Bound<'_, pyo3::types::PyTuple> = key
                .cast()
                .map_err(|_| PyTypeError::new_err("exponent keys must be tuples of ints"))?;
            let exp: Vec<usize> = exp_tuple
                .iter()
                .map(|x| x.extract::<usize>())
                .collect::<PyResult<Vec<_>>>()?;
            if exp.len() != n_vars {
                return Err(PyValueError::new_err(format!(
                    "exponent tuple length {} does not match n_vars={}",
                    exp.len(),
                    n_vars
                )));
            }

            let coeff: Rational = if let Ok(i) = val.extract::<i64>() {
                Rational::new(i, 1)
            } else if let Ok((n, d)) = val.extract::<(i64, i64)>() {
                Rational::new(n, d)
            } else if let Ok(f) = val.extract::<f64>() {
                // Approximate: convert float to rational via continued fractions.
                let bits = 52u32; // f64 mantissa bits
                let scaled = (f * (1i64 << bits) as f64).round() as i64;
                Rational::new(scaled, 1i64 << bits)
            } else {
                return Err(PyTypeError::new_err("coefficients must be int or float"));
            };

            if !domain.is_zero(&coeff) {
                poly_terms.push((exp, coeff));
            }
        }

        Ok(Self {
            inner: SparseMultivariatePolynomial::from_terms(domain, n_vars, poly_terms),
        })
    }

    fn __repr__(&self) -> String {
        format!("MultivariatePolynomial(n_vars={})", self.inner.n_vars())
    }

    fn __str__(&self) -> String {
        format!("MultivariatePolynomial(n_vars={})", self.inner.n_vars())
    }

    fn n_vars(&self) -> usize {
        self.inner.n_vars()
    }
}

/// Convert a list of Python polynomial objects to Rust multivariate polynomials.
///
/// Accepts both `PyPolynomial` (univariate, mapped to variable 0) and
/// `PyMultivariatePolynomial` (native multivariate). The `n_vars` parameter
/// is only used for `PyPolynomial` items; `PyMultivariatePolynomial` carries
/// its own `n_vars`.
fn extract_multivariate_polys(
    polys: &Bound<'_, PyAny>,
    n_vars: usize,
) -> PyResult<Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>> {
    let list: Vec<Bound<'_, PyAny>> = polys
        .extract()
        .map_err(|_| PyTypeError::new_err("expected a list of polynomials"))?;

    if list.is_empty() {
        return Ok(vec![]);
    }

    let mut result = Vec::with_capacity(list.len());
    for item in &list {
        // Try MultivariatePolynomial first (native multivariate).
        if let Ok(mv_poly) = item.extract::<PyRef<'_, PyMultivariatePolynomial>>() {
            result.push(mv_poly.inner.clone());
            continue;
        }
        // Fall back to univariate Polynomial (mapped to variable 0).
        if let Ok(py_poly) = item.extract::<PyRef<'_, PyPolynomial>>() {
            let poly = &py_poly.inner;
            match poly {
                crate::polynomial::PolyErased::Rat(p) => {
                    let terms: Vec<(Vec<usize>, Rational)> = p
                        .coeffs()
                        .iter()
                        .enumerate()
                        .filter(|(_, c)| !RationalDomain.is_zero(c))
                        .map(|(i, c)| {
                            let mut exp = vec![0usize; n_vars];
                            exp[0] = i;
                            (exp, c.clone())
                        })
                        .collect();
                    result.push(SparseMultivariatePolynomial::from_terms(
                        RationalDomain,
                        n_vars,
                        terms,
                    ));
                }
                _ => {
                    return Err(PyTypeError::new_err(
                        "Gröbner basis operations require rational polynomials",
                    ));
                }
            }
        } else {
            return Err(PyTypeError::new_err(
                "expected Polynomial or MultivariatePolynomial objects in the generators list",
            ));
        }
    }
    Ok(result)
}

/// A Gröbner basis computation result.
#[pyclass(name = "GroebnerBasis")]
pub struct PyGroebnerBasis {
    #[pyo3(get)]
    pub n_vars: usize,
    pub basis: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
}

#[pymethods]
impl PyGroebnerBasis {
    fn __len__(&self) -> usize {
        self.basis.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "GroebnerBasis({} elements, {} vars)",
            self.basis.len(),
            self.n_vars
        )
    }

    fn is_groebner_basis(&self) -> bool {
        GroebnerBasis {
            basis: self.basis.clone(),
        }
        .is_groebner_basis()
    }
}

/// A real solution to a polynomial system.
#[pyclass(name = "RealSolution")]
pub struct PyRealSolution {
    #[pyo3(get)]
    pub values: Vec<f64>,
    #[pyo3(get)]
    pub multiplicity: usize,
}

#[pymethods]
impl PyRealSolution {
    fn __repr__(&self) -> String {
        let vals: Vec<String> = self.values.iter().map(|v| format!("{:.6}", v)).collect();
        format!(
            "RealSolution([{}], mult={})",
            vals.join(", "),
            self.multiplicity
        )
    }
}

/// Result of solving a polynomial system.
#[pyclass(name = "PolynomialSystemSolution")]
pub struct PyPolynomialSystemSolution {
    pub inner: PolynomialSystemSolution,
}

#[pymethods]
impl PyPolynomialSystemSolution {
    #[getter]
    fn kind(&self) -> &str {
        match &self.inner {
            PolynomialSystemSolution::ZeroDimensional(_) => "zero_dimensional",
            PolynomialSystemSolution::PositiveDimensional(_) => "positive_dimensional",
            PolynomialSystemSolution::Empty => "empty",
        }
    }

    fn solutions(&self) -> Vec<PyRealSolution> {
        match &self.inner {
            PolynomialSystemSolution::ZeroDimensional(z) => z
                .solutions
                .iter()
                .map(|s| PyRealSolution {
                    values: s.values.clone(),
                    multiplicity: s.multiplicity,
                })
                .collect(),
            _ => vec![],
        }
    }

    #[getter]
    fn vector_space_dimension(&self) -> Option<usize> {
        match &self.inner {
            PolynomialSystemSolution::ZeroDimensional(z) => Some(z.vector_space_dimension),
            _ => None,
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            PolynomialSystemSolution::ZeroDimensional(z) => {
                format!("Solution(zero_dim, {} solutions)", z.solutions.len())
            }
            PolynomialSystemSolution::PositiveDimensional(_) => {
                "Solution(positive_dimensional)".to_string()
            }
            PolynomialSystemSolution::Empty => "Solution(empty)".to_string(),
        }
    }
}

/// A Hilbert series result.
#[pyclass(name = "HilbertSeries")]
pub struct PyHilbertSeries {
    inner: ocas_poly::groebner::hilbert::HilbertSeries,
}

#[pymethods]
impl PyHilbertSeries {
    fn hilbert_function(&self, degree: usize) -> i64 {
        self.inner.hilbert_function(degree)
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    #[getter]
    fn degree(&self) -> i64 {
        self.inner.degree()
    }

    #[getter]
    fn numerator(&self) -> Vec<i64> {
        self.inner.numerator.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "HilbertSeries(dim={}, degree={}, numerator={:?})",
            self.inner.dimension(),
            self.inner.degree(),
            self.inner.numerator
        )
    }
}

/// A primary decomposition component.
#[pyclass(name = "PrimaryComponent")]
pub struct PyPrimaryComponent {
    #[pyo3(get)]
    pub n_vars: usize,
    pub primary: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
    pub prime: Vec<SparseMultivariatePolynomial<RationalDomain, Lex>>,
}

#[pymethods]
impl PyPrimaryComponent {
    fn __repr__(&self) -> String {
        format!(
            "PrimaryComponent(primary={} gens, prime={} gens)",
            self.primary.len(),
            self.prime.len()
        )
    }
}

// ------------------------------------------------------------------
//  Module functions
// ------------------------------------------------------------------

#[pyfunction]
#[pyo3(signature = (generators, n_vars=1, algorithm="auto"))]
pub fn py_groebner_basis(
    generators: &Bound<'_, PyAny>,
    n_vars: usize,
    algorithm: &str,
) -> PyResult<PyGroebnerBasis> {
    let polys = extract_multivariate_polys(generators, n_vars)?;
    let algo = parse_algorithm(algorithm)?;
    let gb = groebner_basis(&polys, algo);
    Ok(PyGroebnerBasis {
        n_vars: polys.first().map(|p| p.n_vars()).unwrap_or(0),
        basis: gb.basis,
    })
}

#[pyfunction]
#[pyo3(signature = (generators, f, n_vars=1, algorithm="auto"))]
pub fn py_ideal_contains(
    generators: &Bound<'_, PyAny>,
    f: &Bound<'_, PyAny>,
    n_vars: usize,
    algorithm: &str,
) -> PyResult<bool> {
    let gens = extract_multivariate_polys(generators, n_vars)?;
    let fs = extract_multivariate_polys(f, n_vars)?;
    let poly = fs
        .into_iter()
        .next()
        .ok_or_else(|| PyValueError::new_err("f must be a non-empty polynomial"))?;
    let algo = parse_algorithm(algorithm)?;
    Ok(ideal::ideal_contains(&gens, &poly, algo))
}

#[pyfunction]
#[pyo3(signature = (equations, n_vars=1, algorithm="auto"))]
pub fn py_solve_polynomial_system(
    equations: &Bound<'_, PyAny>,
    n_vars: usize,
    algorithm: &str,
) -> PyResult<PyPolynomialSystemSolution> {
    let polys = extract_multivariate_polys(equations, n_vars)?;
    let algo = parse_algorithm(algorithm)?;
    let sol = ideal::solve_polynomial_system(&polys, algo);
    Ok(PyPolynomialSystemSolution { inner: sol })
}

#[pyfunction]
pub fn py_hilbert_series(gb: &PyGroebnerBasis) -> PyResult<PyHilbertSeries> {
    let gb_struct = GroebnerBasis {
        basis: gb.basis.clone(),
    };
    let hs = ocas_poly::groebner::hilbert::hilbert_series(&gb_struct);
    Ok(PyHilbertSeries { inner: hs })
}

#[pyfunction]
#[pyo3(signature = (generators, n_vars=1))]
pub fn py_ideal_radical(generators: &Bound<'_, PyAny>, n_vars: usize) -> PyResult<PyGroebnerBasis> {
    let gens = extract_multivariate_polys(generators, n_vars)?;
    let rad = ideal::ideal_radical(&gens);
    Ok(PyGroebnerBasis {
        n_vars: gens.first().map(|p| p.n_vars()).unwrap_or(0),
        basis: rad.basis,
    })
}

#[pyfunction]
#[pyo3(signature = (generators, n_vars=1))]
pub fn py_primary_decomposition(
    generators: &Bound<'_, PyAny>,
    n_vars: usize,
) -> PyResult<Vec<PyPrimaryComponent>> {
    let gens = extract_multivariate_polys(generators, n_vars)?;
    let decomp = ideal::primary_decomposition(&gens);
    let n_vars = gens.first().map(|p| p.n_vars()).unwrap_or(0);
    Ok(decomp
        .into_iter()
        .map(|comp| PyPrimaryComponent {
            n_vars,
            primary: comp.primary,
            prime: comp.prime,
        })
        .collect())
}

#[pyfunction]
pub fn py_is_zero_dimensional(gb: &PyGroebnerBasis) -> bool {
    let gb_struct = GroebnerBasis {
        basis: gb.basis.clone(),
    };
    ideal::is_zero_dimensional(&gb_struct)
}

#[pyfunction]
#[pyo3(signature = (generators, elim_vars, n_vars=1, algorithm="auto"))]
pub fn py_eliminate(
    generators: &Bound<'_, PyAny>,
    elim_vars: usize,
    n_vars: usize,
    algorithm: &str,
) -> PyResult<PyGroebnerBasis> {
    let polys = extract_multivariate_polys(generators, n_vars)?;
    let algo = parse_algorithm(algorithm)?;
    let result = eliminate(&polys, elim_vars, algo);
    Ok(PyGroebnerBasis {
        n_vars: polys.first().map(|p| p.n_vars()).unwrap_or(0),
        basis: result.basis,
    })
}

fn parse_algorithm(s: &str) -> PyResult<Algorithm> {
    match s.to_lowercase().as_str() {
        "auto" => Ok(Algorithm::Auto),
        "f4" => Ok(Algorithm::F4),
        "f5" => Ok(Algorithm::F5),
        "buchberger" => Ok(Algorithm::Buchberger),
        _ => Err(PyValueError::new_err(format!(
            "unknown algorithm '{}': expected auto, f4, f5, or buchberger",
            s
        ))),
    }
}
