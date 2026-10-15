"""Convert PyOptik models into materials for the native solver."""

import numpy as np
import PyOptik

from .material import SellmeierMaterial, SellmeierMedium, TabulatedMaterial, TabulatedMedium
from .units import ureg


def from_pyoptik(model, *, medium=False, extrapolation="error", validate=True):
    """Copy a PyOptik model into a native particle material or surrounding medium.

    PyOptik owns loading and validation. This conversion requests strict index
    validation by default; ``validate=False`` skips sign and passivity checks.
    Media must be lossless under either setting. The solver supports linear
    tables and formula types 1, 2, 5 and 6. Range checking defaults to strict;
    ``extrapolation="linear"`` enables table endpoint slopes or formula
    evaluation beyond declared bounds. Formula validation samples eight points
    within declared bounds; it is not an exhaustive check of the curve.
    """
    if not isinstance(model, (PyOptik.TabulatedMaterial, PyOptik.SellmeierMaterial)):
        raise TypeError("Expected a PyOptik SellmeierMaterial or TabulatedMaterial.")
    if getattr(model, "interpolation", "linear") != "linear":
        raise ValueError("Only linear interpolation is supported by the compiled solver.")
    if extrapolation not in {"error", "linear"}:
        raise ValueError("extrapolation must be 'error' or 'linear'.")
    allow_extrapolation = extrapolation == "linear"
    name = getattr(model, "catalog_id", None) or model.filename
    if isinstance(model, PyOptik.TabulatedMaterial):
        wavelengths, n, k = model.to_arrays(validation="raise" if validate else "ignore")
        indices = n + 1j * k
        if medium and np.any(np.abs(indices.imag) > 1e-12):
            raise ValueError("The native solver requires a lossless medium; extinction coefficients must be zero.")
        model_type = TabulatedMedium if medium else TabulatedMaterial
        return model_type(
            name=name, wavelengths=wavelengths,
            refractive_indices=(indices.real if medium else indices).tolist(),
            allow_extrapolation=allow_extrapolation,
        )

    if model.formula_type not in {1, 2, 5, 6}:
        raise ValueError(f"Formula {model.formula_type} is not supported by the native solver; supported formulas are 1, 2, 5 and 6.")
    bounds = model.wavelength_bound
    bounds = bounds if bounds is not None else [] * ureg.meter
    if validate and len(bounds):
        samples = np.linspace(*bounds.to("meter").magnitude, 8) * ureg.meter
        model.compute_refractive_index(wavelength=samples, out_of_range="raise", validation="raise")
    model_type = SellmeierMedium if medium else SellmeierMaterial
    return model_type(
        name=name, coefficients=model.coefficients, formula_type=model.formula_type,
        wavelength_bound=bounds, allow_extrapolation=allow_extrapolation,
    )


__all__ = ["from_pyoptik"]
