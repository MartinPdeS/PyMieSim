"""Practical material registry, loading, and validation helpers.

This module complements the compiled ``PyMieSim.material`` classes.  It uses
the bundled PyOptik material bank for built-in data and provides a small,
dependency-light file format for user-supplied tabulated indices.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping, Sequence

import numpy as np

from .material import (
    SellmeierMaterial,
    SellmeierMedium,
    TabulatedMaterial,
    TabulatedMedium,
)
from .units import ureg

MaterialKind = Literal["material", "medium", "auto"]
ExtrapolationPolicy = Literal["error", "linear"]
InterpolationPolicy = Literal["linear"]


@dataclass(frozen=True)
class MaterialInfo:
    """Metadata describing one entry in the built-in material registry."""

    name: str
    kind: Literal["sellmeier", "tabulated"]
    wavelength_range: Any
    lossy: bool

    @property
    def wavelength_min(self) -> Any:
        """Lower wavelength covered by the material data."""

        return self.wavelength_range[0]

    @property
    def wavelength_max(self) -> Any:
        """Upper wavelength covered by the material data."""

        return self.wavelength_range[1]


def _pyoptik_bank() -> Any:
    from PyOptik import Material

    # PyMieSim's legacy constructors toggle these filters globally. Reset
    # them so registry queries remain deterministic regardless of prior use.
    Material.use_tabulated = True
    Material.use_sellmeier = True
    return Material


def available_materials(kind: MaterialKind = "auto") -> tuple[str, ...]:
    """Return names in the built-in PyOptik registry.

    Parameters
    ----------
    kind:
        ``"sellmeier"``/``"material"`` returns transparent model entries,
        ``"tabulated"`` returns tabulated entries, and ``"auto"`` returns
        both. ``"material"`` and ``"medium"`` are accepted as convenient
        aliases for ``"auto"`` because the same registry can create either
        a particle material or an optical medium.
    """

    bank = _pyoptik_bank()
    if kind == "tabulated":
        names = bank.tabulated
    elif kind == "sellmeier":
        names = bank.sellmeier
    elif kind in {"auto", "material", "medium"}:
        names = bank.all
    else:
        raise ValueError("kind must be 'auto', 'material', 'medium', 'sellmeier', or 'tabulated'.")
    return tuple(names)


def _as_wavelength_quantity(values: Any, unit: str = "nanometer") -> Any:
    """Convert wavelength values to a one-dimensional metre-compatible quantity."""

    if hasattr(values, "to") and hasattr(values, "magnitude"):
        quantity = values.to("meter")
    else:
        try:
            unit_object = getattr(ureg, unit)
        except AttributeError as error:
            raise ValueError(f"Unknown wavelength unit: {unit!r}") from error
        quantity = np.asarray(values, dtype=float) * unit_object
        quantity = quantity.to("meter")

    magnitudes = np.asarray(quantity.magnitude, dtype=float).reshape(-1)
    if magnitudes.size == 0:
        raise ValueError("At least one wavelength is required.")
    if not np.all(np.isfinite(magnitudes)) or np.any(magnitudes <= 0):
        raise ValueError("Wavelengths must be finite and strictly positive.")
    if magnitudes.size > 1 and np.any(np.diff(magnitudes) <= 0):
        raise ValueError("Wavelengths must be strictly increasing and unique.")
    return magnitudes * ureg.meter


def validate_refractive_indices(
    refractive_indices: Iterable[complex | float],
    *,
    passive: bool = True,
    tolerance: float = 1e-12,
) -> np.ndarray:
    """Validate and return finite complex refractive-index samples.

    PyMieSim uses the optical convention where a passive material has
    ``imag(refractive_index) >= 0``. Negative imaginary parts are therefore
    rejected by default as gain media.
    """

    values = np.asarray(list(refractive_indices), dtype=complex).reshape(-1)
    if values.size == 0:
        raise ValueError("At least one refractive-index sample is required.")
    if not np.all(np.isfinite(values.real)) or not np.all(np.isfinite(values.imag)):
        raise ValueError("Refractive-index samples must be finite.")
    if np.any(values.real <= 0):
        raise ValueError("The real part of the refractive index must be positive.")
    if passive and np.any(values.imag < -tolerance):
        raise ValueError(
            "The material contains gain (negative imaginary refractive index); "
            "passive data requires imag(n) >= 0."
        )
    return values


def validate_tabulated_data(
    wavelengths: Any,
    refractive_indices: Iterable[complex | float],
    *,
    wavelength_unit: str = "nanometer",
    passive: bool = True,
) -> tuple[Any, np.ndarray]:
    """Validate tabulated wavelengths and indices before constructing a model."""

    wavelength_values = _as_wavelength_quantity(wavelengths, wavelength_unit)
    index_values = validate_refractive_indices(refractive_indices, passive=passive)
    if len(wavelength_values.magnitude) != len(index_values):
        raise ValueError("Wavelength and refractive-index arrays must have the same length.")
    return wavelength_values, index_values


def material_info(name: str) -> MaterialInfo:
    """Return metadata for a built-in PyOptik material."""

    if name not in available_materials():
        raise ValueError(f"Unknown material {name!r}. Available materials: {', '.join(available_materials())}")
    model = _pyoptik_bank().get(name)
    if name in _pyoptik_bank().tabulated:
        wavelength_range = model.wavelength_bound.to("meter")
        values = np.asarray(model.n_values) + 1j * np.asarray(model.k_values)
        lossy = bool(np.any(np.abs(np.imag(values)) > 1e-15))
        kind = "tabulated"
    else:
        wavelength_range = model.wavelength_bound.to("meter")
        lossy = False
        kind = "sellmeier"
    return MaterialInfo(name, kind, wavelength_range, lossy)


def _validate_policy(interpolation: str, extrapolation: str) -> None:
    if interpolation != "linear":
        raise ValueError("Only linear interpolation is currently supported by the compiled solver.")
    if extrapolation not in {"error", "linear"}:
        raise ValueError("extrapolation must be 'error' or 'linear'.")


def validate_material(
    model: Any,
    *,
    wavelengths: Any = None,
    passive: bool = True,
    sample_count: int = 8,
) -> None:
    """Validate a material model by sampling its supported wavelength range.

    This catches non-finite, non-positive, or gain-like values in analytical
    models as well as invalid tabulated data. By default the model's declared
    range is used; callers can provide explicit unit-aware wavelengths instead.
    """

    if wavelengths is None:
        if hasattr(model, "wavelengths"):
            wavelengths = model.wavelengths
        elif hasattr(model, "wavelength_bound"):
            wavelengths = model.wavelength_bound
        else:
            raise ValueError("A model must declare wavelengths or wavelength_bound for validation.")

    wavelength_values = _as_wavelength_quantity(wavelengths)
    samples = np.asarray(wavelength_values.magnitude, dtype=float).reshape(-1)
    if samples.size == 2 and sample_count > 2:
        samples = np.linspace(samples[0], samples[-1], sample_count)
    probe = model
    # Older compiled extensions compare Sellmeier bounds in their historical
    # internal unit. Reconstructing a permissive probe lets validation inspect
    # the model values without weakening the policy of the returned model.
    if not getattr(model, "allow_extrapolation", False):
        if hasattr(model, "coefficients") and hasattr(model, "formula_type"):
            model_type = SellmeierMedium if type(model).__name__ == "SellmeierMedium" else SellmeierMaterial
            probe = model_type(
                model.name,
                model.coefficients,
                model.formula_type,
                model.wavelength_bound,
                True,
            )
        elif hasattr(model, "wavelengths") and hasattr(model, "refractive_indices"):
            model_type = TabulatedMedium if type(model).__name__ == "TabulatedMedium" else TabulatedMaterial
            probe = model_type(
                model.name,
                model.wavelengths,
                model.refractive_indices,
                True,
            )
    indices = [probe.get_refractive_index(float(value) * ureg.meter) for value in samples]
    validate_refractive_indices(indices, passive=passive)


def load_material(
    name: str,
    *,
    medium: bool = False,
    interpolation: InterpolationPolicy = "linear",
    extrapolation: ExtrapolationPolicy = "error",
    validate: bool = True,
) -> Any:
    """Load a named material or medium from the built-in registry.

    The default policy is strict range checking. Set ``extrapolation="linear"``
    to use the backend's endpoint-slope linear extrapolation.
    """

    _validate_policy(interpolation, extrapolation)
    bank = _pyoptik_bank()
    if name not in available_materials():
        raise ValueError(f"Unknown material {name!r}. Available materials: {', '.join(available_materials())}")

    source = bank.get(name)
    allow_extrapolation = extrapolation == "linear"
    if name in bank.tabulated:
        wavelengths = source.wavelength.to("meter")
        indices = np.asarray(source.n_values, dtype=float) + 1j * np.asarray(source.k_values, dtype=float)
        if medium:
            indices = indices.real
            if validate:
                validate_refractive_indices(indices, passive=True)
            model = TabulatedMedium(name, wavelengths, indices.real.tolist(), allow_extrapolation)
            if validate:
                validate_material(model)
            return model
        if validate:
            validate_refractive_indices(indices)
        model = TabulatedMaterial(name, wavelengths, indices.tolist(), allow_extrapolation)
        if validate:
            validate_material(model)
        return model

    bounds = source.wavelength_bound.to("meter")
    if medium:
        model = SellmeierMedium(name, source.coefficients, source.formula_type, bounds, allow_extrapolation)
    else:
        model = SellmeierMaterial(name, source.coefficients, source.formula_type, bounds, allow_extrapolation)
    if validate:
        validate_material(model, passive=True)
    return model


def _unit_from_column(name: str, default: str | None = None) -> str | None:
    lowered = name.lower().strip()
    if lowered.endswith("_nm") or lowered.endswith("(nm)"):
        return "nanometer"
    if lowered.endswith("_um") or lowered.endswith("_µm") or lowered.endswith("(um)"):
        return "micrometer"
    if lowered.endswith("_m") or lowered.endswith("(m)"):
        return "meter"
    return default


def _records_from_json(payload: Any) -> tuple[list[Any], list[Any], list[Any]]:
    if isinstance(payload, Mapping):
        records = payload.get("data", payload.get("samples"))
        if records is not None:
            payload = records
        else:
            wavelength = payload.get("wavelength", payload.get("wavelengths"))
            n_values = payload.get("n", payload.get("n_values"))
            if wavelength is None or n_values is None:
                raise ValueError("JSON material data requires wavelength(s) and n/n_values.")
            k_values = payload.get("k", payload.get("k_values", [0.0] * len(wavelength)))
            return list(wavelength), list(n_values), list(k_values)
    if not isinstance(payload, list) or not payload or not isinstance(payload[0], Mapping):
        raise ValueError("JSON material data must be an object or a list of sample objects.")
    wavelength = [row.get("wavelength", row.get("wavelength_nm", row.get("wavelength_um"))) for row in payload]
    n_values = [row.get("n", row.get("real")) for row in payload]
    k_values = [row.get("k", row.get("imag", 0.0)) for row in payload]
    return wavelength, n_values, k_values


def _records_from_csv(path: Path) -> tuple[list[Any], list[Any], list[Any], str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("CSV material data is empty.")
    fields = {field.lower().strip(): field for field in rows[0]}
    wavelength_key = next((fields[key] for key in fields if key in {"wavelength", "wavelength_nm", "wavelength_um", "wavelength_m"}), None)
    n_key = next((fields[key] for key in fields if key in {"n", "real", "refractive_index"}), None)
    k_key = next((fields[key] for key in fields if key in {"k", "imag", "extinction_coefficient"}), None)
    if wavelength_key is None or n_key is None:
        raise ValueError("CSV material data requires wavelength and n columns.")
    unit = _unit_from_column(wavelength_key, "nanometer")
    wavelengths = [float(row[wavelength_key]) for row in rows]
    n_values = [float(row[n_key]) for row in rows]
    k_values = [float(row[k_key]) if k_key and row[k_key] not in {"", None} else 0.0 for row in rows]
    return wavelengths, n_values, k_values, unit


def load_tabulated(
    path: str | Path,
    *,
    name: str | None = None,
    medium: bool = False,
    wavelength_unit: str = "nanometer",
    interpolation: InterpolationPolicy = "linear",
    extrapolation: ExtrapolationPolicy = "error",
    validate: bool = True,
) -> Any:
    """Load a tabulated material from CSV or JSON.

    CSV columns may be ``wavelength,n,k`` or use wavelength suffixes such as
    ``wavelength_nm`` and ``wavelength_um``. JSON accepts either a list of
    ``{"wavelength": ..., "n": ..., "k": ...}`` records or an object with
    array fields ``wavelength(s)``, ``n/n_values``, and optional ``k/k_values``.
    """

    _validate_policy(interpolation, extrapolation)
    source_path = Path(path)
    if source_path.suffix.lower() == ".csv":
        wavelengths, n_values, k_values, detected_unit = _records_from_csv(source_path)
        if detected_unit is not None:
            wavelength_unit = detected_unit
    elif source_path.suffix.lower() == ".json":
        with source_path.open(encoding="utf-8") as stream:
            wavelengths, n_values, k_values = _records_from_json(json.load(stream))
    else:
        raise ValueError("Material data files must have a .csv or .json extension.")

    wavelength_values, index_values = validate_tabulated_data(
        wavelengths,
        np.asarray(n_values, dtype=float) + 1j * np.asarray(k_values, dtype=float),
        wavelength_unit=wavelength_unit,
        passive=validate,
    )
    material_name = name or source_path.stem
    allow_extrapolation = extrapolation == "linear"
    if medium:
        model = TabulatedMedium(material_name, wavelength_values, index_values.real.tolist(), allow_extrapolation)
    else:
        model = TabulatedMaterial(material_name, wavelength_values, index_values.tolist(), allow_extrapolation)
    if validate:
        validate_material(model)
    return model


def validate_wavelength(wavelength: Any, model: Any) -> None:
    """Raise a clear error if a query wavelength is outside a model's range."""

    value = float(wavelength.to("meter").magnitude) if hasattr(wavelength, "to") else float(wavelength)
    if hasattr(model, "wavelengths"):
        bounds = np.asarray(model.wavelengths.to("meter").magnitude, dtype=float)
    elif hasattr(model, "wavelength_bound"):
        bounds = np.asarray(model.wavelength_bound.to("meter").magnitude, dtype=float)
    else:
        return
    if value < bounds[0] or value > bounds[-1]:
        if not getattr(model, "allow_extrapolation", False):
            raise ValueError(
                f"Wavelength {value:g} m is outside the model range "
                f"[{bounds[0]:g}, {bounds[-1]:g}] m and extrapolation is disabled."
            )


__all__ = [
    "ExtrapolationPolicy",
    "InterpolationPolicy",
    "MaterialInfo",
    "MaterialKind",
    "available_materials",
    "load_material",
    "load_tabulated",
    "material_info",
    "validate_refractive_indices",
    "validate_material",
    "validate_tabulated_data",
    "validate_wavelength",
]
