#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyOptik import Material
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import ureg
from PyMieSim.experiment.detector import CoherentModeSet
from PyMieSim.experiment.scatterer import SphereSet
from PyMieSim.experiment.source import GaussianSet, PlaneWaveSet, PolarizationSet
from PyMieSim.experiment import Setup


properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]

measures = SphereSet.available_measure_list


polarization = PolarizationSet(angles=0 * ureg.degree)

gaussian_source = GaussianSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization,
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

planewave_source = PlaneWaveSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization,
    amplitude=[1] * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize(
    "medium_value",
    [[m] for m in medium_properties],
    ids=[f"Medium:{m}" for m in medium_properties],
)

@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources],
)

@pytest.mark.parametrize(
    "refractive_index_value",
    [[m] for m in properties],
    ids=[f"Property:{m}" for m in properties],
)

@pytest.mark.parametrize(
    "measure",
    measures,
)

def test_get_measure(source, measure, refractive_index_value, medium_value):

    kwargs = {}
    for key, value in zip(["", "medium_"], [refractive_index_value, medium_value]):
        if isinstance(value[0], BaseMaterial):
            kwargs[f"{key}material"] = value
        else:
            kwargs[f"{key}refractive_index"] = [value[0].magnitude] * value[0].units

    scatterer = SphereSet(
        diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
        source=source,
        **kwargs
    )

    detector = CoherentModeSet(
        mode_number="LP01",
        rotation=0 * ureg.degree,
        numerical_aperture=[0.1] * ureg.AU,
        polarization_filter=None,
        gamma_offset=[0, 1] * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100,
    )

    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector,
    )

    experiment.get(measure, drop_unique_level=True, scale_unit=True)
    experiment.get(measure, drop_unique_level=False, scale_unit=True, as_numpy=True)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])