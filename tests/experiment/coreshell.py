#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyOptik import Material
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import ureg
from PyMieSim.experiment.detector import PhotodiodeSet
from PyMieSim.experiment.scatterer import CoreShellSet
from PyMieSim.experiment.source import GaussianSet, PlaneWaveSet, PolarizationSet
from PyMieSim.experiment import Setup


core_properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
shell_properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]


polarization_set = PolarizationSet(
    angles=[0] * ureg.radian
)


gaussian_source = GaussianSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2] * ureg.AU,
)

planewave_source = PlaneWaveSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization_set,
    amplitude=[1] * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize(
    "medium_value",
    [[m] for m in medium_properties],
    ids=[f"Medium:{m}" for m in medium_properties]
)

@pytest.mark.parametrize(
    "core_value",
    [[m] for m in core_properties],
    ids=[f"Core:{m}" for m in core_properties]
)

@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources]
)

@pytest.mark.parametrize(
    "shell_value",
    [[m] for m in shell_properties],
    ids=[f"Shell:{m}" for m in shell_properties]
)

@pytest.mark.parametrize(
    "measure",
    CoreShellSet.available_measure_list
)

def test_measure(measure, source, core_value, shell_value, medium_value):
    kwargs = {}

    for key, value in zip(["core", "shell", "medium"], [core_value, shell_value, medium_value]):
        if isinstance(value[0], BaseMaterial):
            kwargs[f"{key}_material"] = value
        else:
            kwargs[f"{key}_refractive_index"] = [value[0].magnitude] * value[0].units

    scatterer = CoreShellSet(
        core_diameter=np.linspace(800, 1000, 10) * ureg.nanometer,
        shell_thickness=[300] * ureg.nanometer,
        source=source,
        **kwargs
    )

    detector = PhotodiodeSet(
        numerical_aperture=[0.2] * ureg.AU,
        gamma_offset=[0] * ureg.degree,
        phi_offset=[0] * ureg.degree,
        sampling=[100],
    )

    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector
    )

    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])