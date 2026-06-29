use pyo3::prelude::*;

/// Python bindings for oCAS.
#[pymodule]
fn ocas_py(m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
