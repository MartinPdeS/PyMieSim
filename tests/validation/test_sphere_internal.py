#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from PyMieSim.single.scatterer import Sphere, CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.mesh import FibonacciMesh  # noqa: F401


def test_Qsca_cross_section():
    """
    Test the consistency between the scattering cross-section obtained directly
    from the sphere object and the one calculated from Qsca and sphere area.
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=1e-6,  # Wavelength in meters
        polarization=0,   # Polarization angle
        optical_power=1,  # Optical power in watts
        NA=0.3            # Numerical aperture
    )

    # Define a spherical scatterer
    sphere = Sphere(
        diameter=300e-9,  # Diameter in meters (e.g., 300 nm)
        index=1.4,        # Refractive index of the sphere
        medium_index=1.0, # Medium index (e.g., air)
        source=source     # Associated light source
    )

    # Calculate scattering cross-section using two different methods
    val0 = sphere.get_cross_section()
    val1 = sphere.Qsca * sphere.area

    # Check if the results are consistent
    assert np.isclose(val0, val1, atol=0, rtol=1e-5), 'Mismatch between cross-section values.'


def test_energy_flow_coupling():
    """
    Test the consistency between the energy flow and coupling values obtained
    from a photodiode detector with a spherical scatterer.
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=1e-6,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    # Define a spherical scatterer
    sphere = Sphere(
        diameter=300e-9,
        index=1.4,
        medium_index=1.0,
        source=source
    )

    # Define a photodiode detector
    detector = Photodiode(
        sampling=500,    # Sampling points for the detector
        NA=2.0,          # Numerical aperture
        gamma_offset=0,  # Offset in the gamma angle
        phi_offset=0     # Offset in the phi angle
    )

    # Calculate energy flow and coupling values
    val0 = detector.get_energy_flow(sphere)
    val1 = detector.coupling(sphere)

    # Check if the results are consistent
    assert np.isclose(val0, val1, atol=0, rtol=1e-5), 'Mismatch between energy flow and coupling values.'


def test_compare_sphere_coreshell_0():
    """
    Compare scattering parameters between a solid sphere and a CoreShell object
    with a zero-thickness shell to verify consistency.
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=1e-6,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    # Define a solid sphere
    sphere = Sphere(
        diameter=1e-6,
        index=1.5,
        source=source,
        medium_index=1.0
    )

    # Define a core-shell scatterer with zero shell thickness
    coreshell = CoreShell(
        core_diameter=1e-6,
        shell_width=0,  # Zero shell width
        core_index=1.5,
        shell_index=1.8,
        medium_index=1.0,
        source=source
    )

    # Compare the scattering parameters between the sphere and core-shell
    for parameter in ['Qsca', 'Qext', 'Qabs']:
        value_sphere = getattr(sphere, parameter)
        value_coreshell = getattr(coreshell, parameter)

        # Check if the results are consistent
        assert np.isclose(value_sphere, value_coreshell, atol=1e-12, rtol=1e-5), (
            f'Mismatch between CoreShell with zero shell and Sphere for parameter: {parameter}'
        )


if __name__ == "__main__":
    pytest.main([__file__])
