#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere, CoreShell
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import IntegratingSphere
from PyMieSim.single import Setup

def test_Qsca_cross_section():
    """
    Test the consistency between the scattering cross-section obtained directly
    from the sphere object and the one calculated from Qsca and sphere area.
    """
    # Define a spherical scatterer
    source = Gaussian(
        wavelength=1000 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    sphere = Sphere(
        diameter=300 * ureg.nanometer,
        material=1.4 * ureg.RIU,
        medium=1.0 * ureg.RIU,
    )

    setup = Setup(
        scatterer=sphere,
        source=source
    )

    # Calculate scattering cross-section using two different methods
    val0 = setup.get("Csca")
    val1 = setup.get("Qsca") * setup.get("cross_section")

    # Check if the results are consistent
    assert np.isclose(
        val0, val1, atol=0, rtol=1e-5
    ), "Mismatch between cross-section values."


def test_energy_flow_coupling():
    """
    Test the consistency between the energy flow and coupling values obtained
    from a photodiode detector with a spherical scatterer.
    """
    source = Gaussian(
        wavelength=1000 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    sphere = Sphere(
        diameter=300 * ureg.nanometer,
        material=1.4 * ureg.RIU,
        medium=1.0 * ureg.RIU,
    )

    detector = IntegratingSphere(
        sampling=500,
    )

    # Calculate energy flow and coupling values
    val0 = detector.get_energy_flow(sphere, source)
    val1 = detector.get_coupling(sphere, source)

    # Check if the results are consistent
    assert np.isclose(
        val0, val1, atol=0, rtol=1e-2
    ), "Mismatch between energy flow and coupling values."


# @pytest.mark.parametrize("parameter", ["Qsca", "Qext", "Qabs"])
# def test_compare_sphere_coreshell_0(parameter: str):
#     """
#     Compare scattering parameters between a solid sphere and a CoreShell object
#     with a zero-thickness shell to verify consistency.
#     """
#     source = Gaussian(
#         wavelength=1000 * ureg.nanometer,
#         polarization=PolarizationState(angle=0 * ureg.degree),
#         optical_power=1 * ureg.watt,
#         numerical_aperture=0.3 * ureg.AU,
#     )

#     sphere = Sphere(
#         diameter=1000 * ureg.nanometer,
#         material=1.5 * ureg.RIU,
#         medium=1.2 * ureg.RIU,
#     )

#     # Define a core-shell scatterer with zero shell thickness
#     coreshell = CoreShell(
#         core_diameter=1000 * ureg.nanometer,
#         shell_thickness=0 * ureg.nanometer,  # Zero shell width
#         core_material=1.5 * ureg.RIU,
#         shell_material=1.8 * ureg.RIU,
#         medium=1.2 * ureg.RIU,
#     )

#     setup_0 = Setup(
#         scatterer=sphere,
#         source=source
#     )

#     setup_1 = Setup(
#         scatterer=coreshell,
#         source=source
#     )

#     # Compare the scattering parameters between the sphere and core-shell
#     value_sphere = setup_0.get(parameter)
#     value_coreshell = setup_1.get(parameter)

#     # Check if the results are consistent
#     assert np.isclose(
#         value_sphere, value_coreshell, atol=1e-12, rtol=1e-5
#     ), f"Mismatch between CoreShell with zero shell and Sphere for parameter: {parameter}"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
