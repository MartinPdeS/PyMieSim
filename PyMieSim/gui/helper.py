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

from typing import Union, Dict, Any
from unittest import case
import numpy
from PyMieSim.experiment.source.gaussian import Gaussian
from PyMieSim.experiment.scatterer.sphere import Sphere
from PyMieSim.experiment.scatterer.cylinder import Cylinder
from PyMieSim.experiment.scatterer.core_shell import CoreShell
from PyMieSim.experiment.detector.photodiode import Photodiode
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, milliwatt, AU, RIU

# Default unit definitions for consistent parameter handling
length_units = nanometer  # Nanometers for all length measurements
power_units = milliwatt   # Milliwatts for optical power measurements
angle_units = degree      # Degrees for angular measurements


def parse_string_to_array_or_float(input_str: str) -> Union[numpy.ndarray, float]:
    """
    Parse a string input to return either a numpy array or a float value.

    This function provides flexible parameter input parsing for the GUI interface,
    supporting multiple input formats to accommodate different use cases in
    Mie scattering experiments.

    Parameters
    ----------
    input_str : str
        Input string in one of the supported formats:
        - 'start:end:count' for linearly spaced arrays (e.g., '400:700:50')
        - 'value1,value2,value3' for comma-separated arrays (e.g., '1.0,1.5,2.0')
        - 'value' for single numeric values (e.g., '1550')

    Returns
    -------
    Union[numpy.ndarray, float]
        - numpy.ndarray: For colon-separated or comma-separated inputs
        - float: For single numeric value inputs

    Raises
    ------
    ValueError
        If the input string format is not recognized or contains invalid
        numeric values that cannot be converted.

    Notes
    -----
    Input Format Examples:
    - Wavelength sweep: '400:800:100' → 100 points from 400 to 800
    - Discrete values: '532,633,1064' → array([532, 633, 1064])
    - Single value: '1550' → 1550.0

    The colon format uses numpy.linspace() for uniform spacing, while
    comma format creates arrays from explicit values. This flexibility
    supports both parametric sweeps and discrete measurement points.

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
    >>> # Single parameter value
    >>> refractive_index = parse_string_to_array_or_float('1.33')
    >>> print(refractive_index)  # 1.33
    """
    try:
        if ":" in input_str:
            # Format: start:end:count → numpy.linspace(start, end, count)
            parts = input_str.split(":")
            if len(parts) != 3:
                raise ValueError("Colon format requires exactly 3 values: 'start:end:count'")
            start, end, count = map(float, parts)
            if count <= 0:
                raise ValueError("Count must be positive")
            return numpy.linspace(start, end, int(count))

        elif "," in input_str:
            # Format: value1,value2,value3 → numpy.array([value1, value2, value3])
            values = [float(val.strip()) for val in input_str.split(",")]
            return numpy.array(values)

        else:
            # Format: single value → float
            return float(input_str.strip())

    except (ValueError, TypeError) as e:
        raise ValueError(
            f"Invalid input string format: '{input_str}'. "
            f"Expected 'start:end:count', comma-separated values, or single numeric value. "
            f"Original error: {str(e)}"
        )


def get_data(
    source_kwargs: Dict[str, str],
    scatterer_kwargs: Dict[str, str],
    detector_kwargs: Dict[str, str],
    measure: str,
    **kwargs: Any
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
    print(scatterer_kwargs)
    try:
        # Create Gaussian beam source with parsed parameters
        source = Gaussian(
            wavelength=parse_string_to_array_or_float(source_kwargs['wavelength']) * length_units,
            polarization=parse_string_to_array_or_float(source_kwargs['polarization']) * angle_units,
            NA=parse_string_to_array_or_float(source_kwargs['NA']) * AU,
            optical_power=parse_string_to_array_or_float(source_kwargs['optical_power']) * power_units
        )

        # Create spherical scatterer with parsed parameters

        match scatterer_kwargs['type']:
            case 'sphere':
                scatterer = Sphere(
                    diameter=parse_string_to_array_or_float(scatterer_kwargs['diameter']) * length_units,
                    property=parse_string_to_array_or_float(scatterer_kwargs['property']) * RIU,
                    medium_property=parse_string_to_array_or_float(scatterer_kwargs['medium_property']) * RIU,
                    source=source
                )
            case 'coreshell':
                scatterer = CoreShell(
                    core_diameter=parse_string_to_array_or_float(scatterer_kwargs['core_diameter']) * length_units,
                    shell_thickness=parse_string_to_array_or_float(scatterer_kwargs['shell_thickness']) * length_units,
                    core_property=parse_string_to_array_or_float(scatterer_kwargs['core_property']) * RIU,
                    shell_property=parse_string_to_array_or_float(scatterer_kwargs['shell_property']) * RIU,
                    medium_property=parse_string_to_array_or_float(scatterer_kwargs['medium_property']) * RIU,
                    source=source
                )
            case 'cylinder':
                scatterer = Cylinder(
                    diameter=parse_string_to_array_or_float(scatterer_kwargs['diameter']) * length_units,
                    property=parse_string_to_array_or_float(scatterer_kwargs['property']) * RIU,
                    medium_property=parse_string_to_array_or_float(scatterer_kwargs['medium_property']) * RIU,
                    source=source
                )

        # Create photodiode detector with parsed parameters
        detector = Photodiode(
            NA=parse_string_to_array_or_float(detector_kwargs['NA']) * AU,
            gamma_offset=parse_string_to_array_or_float(detector_kwargs['gamma_offset']) * angle_units,
            phi_offset=parse_string_to_array_or_float(detector_kwargs['phi_offset']) * angle_units,
            sampling=parse_string_to_array_or_float(detector_kwargs['sampling']) * AU,
        )

        # Set up the complete experiment
        setup = Setup(source=source, scatterer=scatterer, detector=detector)

        # Execute the measurement and return results
        dataframe = setup.get(measure, **kwargs)

        return dataframe

    except KeyError as e:
        raise KeyError(f"Missing required parameter: {str(e)}. Check that all required keys are present in parameter dictionaries.")
    except Exception as e:
        raise RuntimeError(f"Failed to execute Mie scattering experiment: {str(e)}. Check parameter values and measurement type.")
