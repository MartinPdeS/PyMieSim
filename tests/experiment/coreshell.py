#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyOptik import Material
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import ureg
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import CoreShell
from PyMieSim.experiment.source import Gaussian, PlaneWave, PolarizationSet
from PyMieSim.experiment import Setup


core_properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
shell_properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]


polarization_set = PolarizationSet(
    angles=[0] * ureg.radian
)


gaussian_source = Gaussian(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

planewave_source = PlaneWave(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=polarization_set,
    amplitude=1 * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


def normalize_material(value):
    """
    Always return a list of BaseMaterial.
    """
    if isinstance(value, (list, tuple, np.ndarray)):
        return list(value)

    return [value]


def normalize_refractive_index(value):
    """
    Always return a Pint quantity with vector magnitude.
    """
    if not hasattr(value, "units"):
        raise TypeError("Refractive index must be a Pint quantity")

    return np.atleast_1d(value.magnitude) * value.units


def build_core_shell_kwargs(core_value, shell_value, medium_value):

    kwargs = {}

    if isinstance(core_value, BaseMaterial):
        kwargs["core_material"] = normalize_material(core_value)
    else:
        kwargs["core_refractive_index"] = normalize_refractive_index(core_value)

    if isinstance(shell_value, BaseMaterial):
        kwargs["shell_material"] = normalize_material(shell_value)
    else:
        kwargs["shell_refractive_index"] = normalize_refractive_index(shell_value)

    if isinstance(medium_value, BaseMaterial):
        kwargs["medium_material"] = normalize_material(medium_value)
    else:
        kwargs["medium_refractive_index"] = normalize_refractive_index(medium_value)

    return kwargs


@pytest.mark.parametrize(
    "medium_value",
    medium_properties,
    ids=[f"Medium:{m}" for m in medium_properties]
)

@pytest.mark.parametrize(
    "core_value",
    core_properties,
    ids=[f"Core:{m}" for m in core_properties]
)

@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources]
)

@pytest.mark.parametrize(
    "shell_value",
    shell_properties,
    ids=[f"Shell:{m}" for m in shell_properties]
)

@pytest.mark.parametrize(
    "measure",
    CoreShell.available_measure_list
)

def test_measure(measure, source, core_value, shell_value, medium_value):

    scatterer_kwargs = build_core_shell_kwargs(
        core_value,
        shell_value,
        medium_value
    )

    print(scatterer_kwargs)

    scatterer = CoreShell(
        core_diameter=np.linspace(800, 1000, 10) * ureg.nanometer,
        shell_thickness=[300] * ureg.nanometer,
        source=source,
        **scatterer_kwargs
    )

    detector = Photodiode(
        numerical_aperture=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100,
    )

    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector
    )

    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])