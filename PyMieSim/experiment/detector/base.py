#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from pint_pandas import PintArray
from pydantic import field_validator
from PyMieSim.binary.interface_sets import CppDetectorSet
from PyMieSim.units import Quantity, radian, degree


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
    NA : Quantity
        The numerical aperture of the detector, typically a float-like value wrapped in a physical unit (e.g., steradian).
    gamma_offset : Quantity
        The angular offset in the gamma direction (in degrees).
    phi_offset : Quantity
        The angular offset in the phi direction (in degrees).
    rotation : Quantity
        The rotational angle of the detector (in degrees).
    mean_coupling : bool
        Specifies whether the mean coupling of detected modes is to be considered (default is False).
    coherent : bool
        Specifies whether the detection process is coherent (default is True).
    sampling : Optional[Quantity]
        The sampling rate of the detector. If not specified, defaults to 200.
    polarization_filter : Optional[Quantity]
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
    @field_validator('mode_number', mode='plain')
    def _validate_mode_number(cls, mode_number):
        """Ensure mode numbers are valid and belong to supported families."""
        mode_number = numpy.atleast_1d(mode_number).astype(str)
        for mode in mode_number:
            if mode[:2] not in ['LP', 'HG', 'LG', 'NC']:
                raise ValueError(f'Invalid mode family {mode[:2]}. Must be one of: LP, HG, LG, NC')
        return mode_number

    @field_validator('polarization_filter', mode='plain')
    def validate_polarization(cls, value):
        """
        Validates the polarization filter. If not provided, defaults to NaN degrees.
        Ensures the value has angular units (degree or radian).

        Parameters
        ----------
        value : Any
            The input value for polarization filter.

        Returns
        -------
        numpy.ndarray
            A NumPy array containing the validated and converted polarization filter value.
        """
        if value is None:
            value = numpy.nan * degree

        if not isinstance(value, Quantity) or not value.check(degree):
            raise ValueError(f"{value} must have angle units (degree or radian).")

        return numpy.atleast_1d(value).astype(float)

    @field_validator('gamma_offset', 'phi_offset', 'rotation', mode='plain')
    def validate_angle_quantity(cls, value):
        """
        Ensures that angular quantities (gamma_offset, phi_offset, rotation) are correctly formatted as Quantities with angle units.

        Parameters
        ----------
        value : Any
            The input value for the angle.

        Returns
        -------
        numpy.ndarray
            A NumPy array containing the validated and converted angle value.
        """
        if not isinstance(value, Quantity) or not value.check(degree):
            raise ValueError(f"{value} must be a Quantity with angle units [degree or radian].")

        return numpy.atleast_1d(value)

    @field_validator('NA', 'cache_NA', 'sampling', mode='plain')
    def validate_au_quantity(cls, value):
        """
        Ensures that numerical values such as numerical aperture (NA) and sampling rate are correctly cast into NumPy arrays.

        Parameters
        ----------
        value : Any
            The input value to be validated.

        Returns
        -------
        numpy.ndarray
            A NumPy array representing the validated input value.
        """
        if not isinstance(value, Quantity) or not value.check(degree):
            raise ValueError(f"{value} must be a Quantity with arbitrary units [AU].")

        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    def _generate_binding(self) -> None:
        """
        Initializes the C++ binding for the detector using the given simulation parameters. This ensures that the
        detector is correctly linked to the backend, enabling high-performance Mie scattering calculations.

        Sets up parameters such as mode number, sampling rate, NA, and various offsets for the simulation.
        """
        self.binding_kwargs = {
            "mode_number": self.mode_number,
            "sampling": self.sampling,
            "NA": self.NA,
            "cache_NA": self.cache_NA,
            "polarization_filter": self.polarization_filter.to(radian).magnitude if self.polarization_filter is not None else numpy.nan,
            "phi_offset": self.phi_offset.to(radian).magnitude,
            "gamma_offset": self.gamma_offset.to(radian).magnitude,
            "rotation": self.rotation.to(radian).magnitude,
            "is_sequential": self.is_sequential
        }

        # Ensure all values are at least 1D arrays for compatibility
        self.binding_kwargs = {k: numpy.atleast_1d(v) for k, v in self.binding_kwargs.items()}

        # Additional detector settings
        self.binding_kwargs["mean_coupling"] = self.mean_coupling
        self.binding_kwargs["coherent"] = self.coherent

        self.binding = CppDetectorSet(**self.binding_kwargs)

    def _generate_mapping(self) -> None:
        """
        Updates the internal mapping of the detector with current parameter values, allowing for visual representation
        of the detector's properties in a tabular format (useful for debugging and visualization).

        Attributes like mode number, NA, offsets, and sampling are included in this mapping.
        """
        self.mapping = {}
        self.mapping.update({'detector:mode_number': self.mode_number})

        for attr in ['NA', 'sampling', 'cache_NA', 'phi_offset', 'gamma_offset', 'rotation', 'polarization_filter']:
            self.mapping["detector:" + attr] = PintArray(getattr(self, attr), dtype=getattr(self, attr).units)
