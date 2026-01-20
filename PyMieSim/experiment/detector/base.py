#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from pint_pandas import PintArray
from pydantic import field_validator
from TypedUnit import ureg, Angle, Dimensionless

class BaseDetector:
    """
    A base class for defining detectors in Mie scattering simulations.

    This class provides the foundational structure for defining detectors that interface with C++ backends to perform
    high-performance scattering calculations. It is designed to handle common parameters such as numerical aperture (NA),
    angular offsets, rotation, and polarization. Subclasses should extend this base class to create specific detector types.

    Attributes
    ----------
    mode_number : str
        The mode number used in the detection process, identifying the angular momentum mode involved.
    NA : Dimensionless
        The numerical aperture of the detector, typically a float-like value wrapped in a physical unit (e.g., steradian).
    gamma_offset : Angle
        The angular offset in the gamma direction (in degrees).
    phi_offset : Angle
        The angular offset in the phi direction (in degrees).
    rotation : Angle
        The rotational angle of the detector (in degrees).
    mean_coupling : bool
        Specifies whether the mean coupling of detected modes is to be considered (default is False).
    coherent : bool
        Specifies whether the detection process is coherent (default is True).
    sampling : Optional[Dimensionless]
        The sampling rate of the detector. If not specified, defaults to 200.
    polarization_filter : Optional[Angle]
        The polarization filter angle (in degrees). Defaults to NaN if not provided.

    Methods
    -------
    validate_polarization_filter(cls, value)
        Ensures that the polarization filter is either None or a Quantity with angle units (degree or radian).
    validate_angle_quantity(cls, value)
        Validates that angles like gamma_offset, phi_offset, and rotation are Quantities with angle units.
    validate_au_quantity(cls, value)
        Ensures that numerical values are correctly cast as NumPy arrays.
    __post_init__()
        Initializes the internal detector mapping and binds the detector to the C++ backend.
    _initialize_binding()
        Sets up the C++ bindings required for efficient simulation by passing detector parameters.
    _generate_mapping()
        Generates a mapping of scatterer properties to be used for visualization purposes.
    """
    def _generate_mapping(self) -> None:
        """
        Updates the internal mapping of the detector with current parameter values, allowing for visual representation
        of the detector's properties in a tabular format (useful for debugging and visualization).

        Attributes like mode number, NA, offsets, and sampling are included in this mapping.
        """
        self.mapping = {}

        for attr in self.attributes:
            value = getattr(self, attr)
            string = "detector:" + attr
            if hasattr(value, "magnitude"):
                self.mapping[string] = PintArray(value, dtype=value.units)
            else:
                self.mapping[string] = [repr(v) for v in value]