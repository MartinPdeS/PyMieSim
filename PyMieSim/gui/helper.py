#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PyMieSim GUI Helper Module

This module provides utility functions for the PyMieSim graphical user interface,
enabling easy setup and execution of Mie scattering experiments through string-based
parameter input and automated data processing.

Key functionalities:
- String-to-array parsing for flexible parameter input
- Automated experiment setup with source, scatterer, and detector configuration
- Data retrieval and processing for various scattering measurements
- Unit conversion and validation for physical parameters

The module serves as a bridge between the GUI interface and the core PyMieSim
computational engine, handling parameter parsing, validation, and experiment
orchestration.
"""

from typing import Union, Dict, Any, List
import numpy
import re
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment.scatterer import Sphere, Cylinder, CoreShell
from PyMieSim.experiment.detector import Photodiode, CoherentMode
from PyMieSim.experiment import Setup
from TypedUnit import ureg

# Default unit definitions for consistent parameter handling
length_units = ureg.nanometer  # Nanometers for all length measurements
power_units = ureg.milliwatt  # Milliwatts for optical power measurements
angle_units = ureg.degree  # Degrees for angular measurements


def parse_string_to_array_or_float(
    input_str: str,
) -> Union[numpy.ndarray, float, List[str], complex]:
    """
    Parse a string input to return either a numpy array, float value, string list, or complex number.

    This function provides flexible parameter input parsing for the GUI interface,
    supporting multiple input formats to accommodate different use cases in
    Mie scattering experiments.

    Parameters
    ----------
    input_str : str
        Input string in one of the supported formats:
        - 'start:end:count' for linearly spaced arrays (e.g., '400:700:50')
        - 'value1,value2,value3' for comma-separated arrays (e.g., '1.0,1.5,2.0')
        - 'str1,str2,str3' for comma-separated string lists (e.g., 'LP01,HG11')
        - 'real+imagj' for complex numbers (e.g., '1.5+0.1j' or '2-0.5j')
        - 'value' for single numeric values (e.g., '1550')

    Returns
    -------
    Union[numpy.ndarray, float, List[str], complex]
        - numpy.ndarray: For colon-separated or comma-separated numeric inputs
        - float: For single numeric value inputs
        - List[str]: For comma-separated string inputs
        - complex: For complex number inputs

    Raises
    ------
    ValueError
        If the input string format is not recognized or contains invalid
        values that cannot be converted.

    Notes
    -----
    Input Format Examples:
    - Wavelength sweep: '400:800:100' → 100 points from 400 to 800
    - Discrete values: '532,633,1064' → array([532, 633, 1064])
    - Complex values: '1.5+0.1j,2.0-0.3j' → array([1.5+0.1j, 2.0-0.3j])
    - String list: 'LP01,HG11,TM01' → ['LP01', 'HG11', 'TM01']
    - Single value: '1550' → 1550.0
    - Complex single: '1.5+0.1j' → (1.5+0.1j)

    The colon format uses numpy.linspace() for uniform spacing, while
    comma format creates arrays from explicit values. This flexibility
    supports both parametric sweeps, discrete measurement points, mode
    labels, and complex refractive indices.

    Examples
    --------
    >>> # Linear sweep for wavelength range
    >>> wavelengths = parse_string_to_array_or_float('400:800:100')
    >>> print(wavelengths.shape)  # (100,)
    >>>
    >>> # Discrete measurement points
    >>> diameters = parse_string_to_array_or_float('100,200,500,1000')
    >>> print(diameters)  # [100. 200. 500. 1000.]
    >>>
    >>> # Mode labels
    >>> modes = parse_string_to_array_or_float('LP01,HG11,TM01')
    >>> print(modes)  # ['LP01', 'HG11', 'TM01']
    >>>
    >>> # Complex refractive index
    >>> n_complex = parse_string_to_array_or_float('1.5+0.1j')
    >>> print(n_complex)  # (1.5+0.1j)
    >>>
    >>> # Complex array
    >>> n_array = parse_string_to_array_or_float('1.5+0.1j,2.0-0.3j')
    >>> print(n_array)  # [1.5+0.1j 2.0-0.3j]
    """

    def _parse_complex_or_float(value_str: str) -> Union[float, complex]:
        """Helper function to parse a single value as complex or float."""
        value_str = value_str.strip()

        # Check if it's a complex number (contains 'j' or 'J')
        if "j" in value_str.lower():
            try:
                # Handle Python complex number format
                return complex(value_str.replace("i", "j").replace("I", "j"))
            except ValueError:
                raise ValueError(f"Invalid complex number format: '{value_str}'")
        else:
            try:
                return float(value_str)
            except ValueError:
                raise ValueError(f"Invalid numeric format: '{value_str}'")

    def _is_numeric_value(value_str: str) -> bool:
        """Check if a string represents a numeric value (float or complex)."""
        value_str = value_str.strip()

        # Check for complex number pattern
        complex_pattern = r"^[+-]?(\d+\.?\d*|\.\d+)([+-]\d+\.?\d*[ji])?$"
        if re.match(complex_pattern, value_str.lower()):
            return True

        # Check for simple float
        try:
            float(value_str)
            return True
        except ValueError:
            return False

    try:
        input_str = input_str.strip()

        if ":" in input_str:
            # Format: start:end:count → numpy.linspace(start, end, count)
            parts = input_str.split(":")
            if len(parts) != 3:
                raise ValueError(
                    "Colon format requires exactly 3 values: 'start:end:count'"
                )
            start, end, count = map(float, parts)
            if count <= 0:
                raise ValueError("Count must be positive")
            return numpy.linspace(start, end, int(count))

        elif "," in input_str:
            # Format: value1,value2,value3 → determine if numeric or string
            values = [val.strip() for val in input_str.split(",")]

            # Check if all values are numeric
            if all(_is_numeric_value(val) for val in values):
                # Parse as numeric values (float or complex)
                parsed_values = [_parse_complex_or_float(val) for val in values]
                return numpy.array(parsed_values)
            else:
                # Return as string list
                return values

        else:
            # Single value: determine if numeric, complex, or string
            if _is_numeric_value(input_str):
                return _parse_complex_or_float(input_str)
            else:
                # Return as single string (could be a mode label like 'LP01')
                return input_str

    except (ValueError, TypeError) as e:
        raise ValueError(
            f"Invalid input string format: '{input_str}'. "
            f"Expected 'start:end:count', comma-separated values (numeric or string), "
            f"complex numbers, or single numeric/string value. "
            f"Original error: {str(e)}"
        )


def get_data(
    source_kwargs: Dict[str, str],
    scatterer_kwargs: Dict[str, str],
    detector_kwargs: Dict[str, str],
    measure: str,
    **kwargs: Any,
) -> Any:
    """
    Execute a complete Mie scattering experiment and retrieve measurement data.

    This function orchestrates a full PyMieSim experiment by setting up the source,
    scatterer, and detector based on string-formatted parameters, then computing
    the requested scattering measurement. It serves as the main interface between
    the GUI and the computational engine.

    Parameters
    ----------
    source_kwargs : Dict[str, str]
        Dictionary containing source parameters as strings:
        - 'wavelength': Incident wavelength(s) in nanometers
        - 'polarization': Polarization angle(s) in degrees (0° = x-polarized)
        - 'NA': Numerical aperture (dimensionless, 0-1)
        - 'optical_power': Incident optical power in milliwatts

    scatterer_kwargs : Dict[str, str]
        Dictionary containing scatterer parameters as strings:
        - 'diameter': Particle diameter(s) in nanometers
        - 'property': Particle refractive index (complex or real)
        - 'medium_property': Medium refractive index (typically 1.0 for air)

    detector_kwargs : Dict[str, str]
        Dictionary containing detector parameters as strings:
        - 'NA': Detector numerical aperture (collection efficiency)
        - 'gamma_offset': Polar angle offset in degrees
        - 'phi_offset': Azimuthal angle offset in degrees
        - 'sampling': Number of angular sampling points

    measure : str
        Type of measurement to perform. Common options include:
        - 'Qsca': Scattering efficiency
        - 'Qext': Extinction efficiency
        - 'Qabs': Absorption efficiency
        - 'Qback': Backscattering efficiency
        - 'g': Asymmetry parameter
        - 'coupling': Coupling efficiency to detector

    **kwargs : Any
        Additional keyword arguments passed to the Setup.get() method
        for advanced measurement configurations.

    Returns
    -------
    Any
        Measurement results as a pandas DataFrame or numpy array,
        depending on the measurement type and parameter dimensions.
        The structure varies based on:
        - Single vs. multiple parameter values
        - Measurement type (efficiency, cross-section, field)
        - Additional processing options in kwargs

    Raises
    ------
    ValueError
        If parameter strings cannot be parsed or contain invalid values.
    KeyError
        If required parameter keys are missing from input dictionaries.
    RuntimeError
        If the computational setup fails or measurement cannot be completed.

    """
    match source_kwargs["type"]:
        case "planewave":
            source = PlaneWave(
                wavelength=parse_string_to_array_or_float(source_kwargs["wavelength"])
                * length_units,
                polarization=parse_string_to_array_or_float(
                    source_kwargs["polarization"]
                )
                * angle_units,
                amplitude=parse_string_to_array_or_float(source_kwargs["amplitude"])
                * ureg.volt
                / ureg.meter,
            )
        case "gaussian":
            source = Gaussian(
                wavelength=parse_string_to_array_or_float(source_kwargs["wavelength"])
                * length_units,
                polarization=parse_string_to_array_or_float(
                    source_kwargs["polarization"]
                )
                * angle_units,
                NA=parse_string_to_array_or_float(source_kwargs["numerical_aperture"])
                * ureg.AU,
                optical_power=parse_string_to_array_or_float(
                    source_kwargs["optical_power"]
                )
                * power_units,
            )

    match scatterer_kwargs["type"]:
        case "sphere":
            scatterer = Sphere(
                diameter=parse_string_to_array_or_float(scatterer_kwargs["diameter"])
                * length_units,
                property=parse_string_to_array_or_float(scatterer_kwargs["property"])
                * ureg.RIU,
                medium_property=parse_string_to_array_or_float(
                    scatterer_kwargs["medium_property"]
                )
                * ureg.RIU,
                source=source,
            )
        case "coreshell":
            scatterer = CoreShell(
                core_diameter=parse_string_to_array_or_float(
                    scatterer_kwargs["core_diameter"]
                )
                * length_units,
                shell_thickness=parse_string_to_array_or_float(
                    scatterer_kwargs["shell_thickness"]
                )
                * length_units,
                core_property=parse_string_to_array_or_float(
                    scatterer_kwargs["core_property"]
                )
                * ureg.RIU,
                shell_property=parse_string_to_array_or_float(
                    scatterer_kwargs["shell_property"]
                )
                * ureg.RIU,
                medium_property=parse_string_to_array_or_float(
                    scatterer_kwargs["medium_property"]
                )
                * ureg.RIU,
                source=source,
            )
        case "cylinder":
            scatterer = Cylinder(
                diameter=parse_string_to_array_or_float(scatterer_kwargs["diameter"])
                * length_units,
                property=parse_string_to_array_or_float(scatterer_kwargs["property"])
                * ureg.RIU,
                medium_property=parse_string_to_array_or_float(
                    scatterer_kwargs["medium_property"]
                )
                * ureg.RIU,
                source=source,
            )

    match detector_kwargs["type"]:
        case "photodiode":
            detector = Photodiode(
                NA=parse_string_to_array_or_float(detector_kwargs["numerical_aperture"])
                * ureg.AU,
                gamma_offset=parse_string_to_array_or_float(
                    detector_kwargs["gamma_offset"]
                )
                * angle_units,
                phi_offset=parse_string_to_array_or_float(detector_kwargs["phi_offset"])
                * angle_units,
                sampling=parse_string_to_array_or_float(detector_kwargs["sampling"])
                * ureg.AU,
            )
        case "coherentmode":
            detector = CoherentMode(
                mode_number=parse_string_to_array_or_float(
                    detector_kwargs["mode_number"]
                ),
                NA=parse_string_to_array_or_float(detector_kwargs["numerical_aperture"])
                * ureg.AU,
                gamma_offset=parse_string_to_array_or_float(
                    detector_kwargs["gamma_offset"]
                )
                * angle_units,
                phi_offset=parse_string_to_array_or_float(detector_kwargs["phi_offset"])
                * angle_units,
                sampling=parse_string_to_array_or_float(detector_kwargs["sampling"])
                * ureg.AU,
                rotation=parse_string_to_array_or_float(
                    detector_kwargs.get("rotation", 0)
                )
                * angle_units,
            )

    # Set up the complete experiment
    setup = Setup(source=source, scatterer=scatterer, detector=detector)

    # Execute the measurement and return results
    dataframe = setup.get(measure, **kwargs)

    return dataframe
