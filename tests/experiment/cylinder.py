#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyOptik import Material
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import ureg
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import InfiniteCylinder
from PyMieSim.experiment.source import Gaussian, PlaneWave, PolarizationSet
from PyMieSim.experiment import Setup


properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]

measures = InfiniteCylinder.available_measure_list


polarization = PolarizationSet(angles=0 * ureg.degree)

gaussian_source = Gaussian(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=polarization,
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

planewave_source = PlaneWave(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=polarization,
    amplitude=1 * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


def normalize_material(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        return list(value)

    return [value]


def normalize_refractive_index(value):
    return np.atleast_1d(value.magnitude) * value.units


def build_kwargs(refractive_index_value, medium_value):

    kwargs = {}

    if isinstance(refractive_index_value, BaseMaterial):
        kwargs["material"] = normalize_material(refractive_index_value)
    else:
        kwargs["refractive_index"] = normalize_refractive_index(refractive_index_value)

    if isinstance(medium_value, BaseMaterial):
        kwargs["medium_material"] = normalize_material(medium_value)
    else:
        kwargs["medium_refractive_index"] = normalize_refractive_index(medium_value)

    return kwargs


@pytest.mark.parametrize(
    "medium_value",
    medium_properties,
    ids=[f"Medium:{m}" for m in medium_properties],
)

@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources],
)

@pytest.mark.parametrize(
    "refractive_index_value",
    properties,
    ids=[f"Property:{m}" for m in properties],
)

@pytest.mark.parametrize(
    "measure",
    measures,
)

def test_measure(measure, source, medium_value, refractive_index_value):

    scatterer_kwargs = build_kwargs(
        refractive_index_value,
        medium_value,
    )

    scatterer = InfiniteCylinder(
        diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
        source=source,
        **scatterer_kwargs
    )

    detector = Photodiode(
        numerical_aperture=0.2 * ureg.AU,
        polarization_filter=None,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100,
    )

    experiment = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector,
    )

    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])