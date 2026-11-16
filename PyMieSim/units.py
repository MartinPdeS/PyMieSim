"""Shared TypedUnit registry, physical types, and validation API.

Validation preserves supplied quantities. Native bindings convert magnitudes
only after TypedUnit has checked units, shape, and finiteness.
"""

from TypedUnit import (
    Angle,
    Dimensionless,
    ElectricField,
    Finite,
    Length,
    Matrix,
    NoUnit,
    Power,
    Quantity,
    RefractiveIndex,
    Scalar,
    Vector,
    unit_validator,
    ureg,
    validate_annotation,
    validate_quantity,
)
from PyMieSim import _pint

_pint.set_ureg(ureg)
