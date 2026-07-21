"""Parsing helpers for dashboard form values."""

from __future__ import annotations

from typing import Any

import numpy as np

from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.material import SellmeierMaterial, SellmeierMedium, TabulatedMaterial


def parse_expression(raw_value: Any) -> Any:
    """Parse a scalar, CSV sequence, or ``start:end:count`` specification."""
    if raw_value is None:
        return None

    if isinstance(raw_value, (int, float, complex, np.ndarray, list, tuple)):
        return raw_value

    text = str(raw_value).strip()

    if text == "":
        return None

    if text.count(":") == 2:
        start_text, stop_text, count_text = [part.strip() for part in text.split(":")]
        count = int(float(count_text))
        if count <= 0:
            raise ValueError("Range count must be positive.")
        return np.linspace(float(start_text), float(stop_text), count)

    if "," in text:
        tokens = [_parse_scalar_or_text(token.strip()) for token in text.split(",") if token.strip()]

        if all(not isinstance(token, str) for token in tokens):
            dtype = complex if any(isinstance(token, complex) and not isinstance(token, bool) for token in tokens) else float
            return np.asarray(tokens, dtype=dtype)

        return [str(token) for token in tokens]

    return _parse_scalar_or_text(text)


def parse_numeric_expression(raw_value: Any, *, integer: bool = False) -> Any:
    """Parse a numeric scalar or numeric sequence."""
    parsed = parse_expression(raw_value)

    if parsed is None:
        return None

    if isinstance(parsed, str):
        raise ValueError(f"Expected numeric values, received '{parsed}'.")

    if isinstance(parsed, np.ndarray):
        return parsed.astype(int if integer else float)

    if isinstance(parsed, list):
        if any(isinstance(value, str) for value in parsed):
            raise ValueError(f"Expected numeric values, received '{parsed}'.")
        return np.asarray(parsed, dtype=int if integer else float)

    return int(parsed) if integer else float(parsed)


def parse_quantity_expression(raw_value: Any, unit: Any) -> Any:
    """Parse a numeric expression and attach a Pint unit."""
    parsed = parse_numeric_expression(raw_value)

    if parsed is None:
        return None

    return parsed * unit


def parse_polarization(raw_value: Any, unit: Any) -> PolarizationSet:
    """Parse a polarization angle expression into a ``PolarizationSet``."""
    angles = parse_quantity_expression(raw_value, unit)
    return PolarizationSet(angles=angles)


def parse_mode_numbers(raw_value: Any) -> Any:
    """Parse one or more coherent mode labels."""
    parsed = parse_expression(raw_value)

    if isinstance(parsed, list):
        return parsed

    return parsed


def parse_material_values(raw_value: Any, *, medium: bool = False) -> Any:
    """Parse numeric or named material or medium definitions."""
    parsed = parse_expression(raw_value)

    if parsed is None:
        return None

    if isinstance(parsed, np.ndarray):
        return parsed

    if isinstance(parsed, list):
        return [_resolve_material_entry(value, medium=medium) for value in parsed]

    return _resolve_material_entry(parsed, medium=medium)


def serialize_value(value: Any) -> Any:
    """Convert NumPy and object values into JSON-friendly primitives."""
    if isinstance(value, np.generic):
        return value.item()

    if isinstance(value, complex):
        return str(value)

    if isinstance(value, (int, float, str, bool)) or value is None:
        return value

    return str(value)


def _parse_scalar_or_text(text: str) -> Any:
    """Parse one scalar token as float, complex, or raw text."""
    if text == "":
        return None

    try:
        if "j" in text.lower():
            return complex(text.replace("i", "j").replace("I", "j"))
        return float(text)
    except ValueError:
        return text


def _resolve_material_entry(value: Any, *, medium: bool) -> Any:
    """Resolve one material token into either a numeric value or material object."""
    if not isinstance(value, str):
        return value

    constructors = (SellmeierMedium,) if medium else (SellmeierMaterial, TabulatedMaterial)
    last_error = None

    for constructor in constructors:
        try:
            return constructor(value)
        except Exception as error:  # pragma: no cover
            last_error = error

    raise ValueError(f"Unknown {'medium' if medium else 'material'} '{value}'.") from last_error
