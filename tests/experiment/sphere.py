#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from PyMieSim.units import ureg
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.material import SellmeierMaterial, TabulatedMaterial, SellmeierMedium


properties = [
    [TabulatedMaterial("silver")],
    # [SellmeierMaterial("fused_silica")],
    # [1.4] * ureg.RIU
]

medium_properties = [
    [SellmeierMedium("water")],
    # [1.1] * ureg.RIU
]

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
    "medium", medium_properties,
    ids=[f"Medium:{m}" for m in medium_properties],
)

@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources],
)

@pytest.mark.parametrize(
    "material", properties,
    ids=[f"Property:{m}" for m in properties],
)

@pytest.mark.parametrize("measure", measures)
def test_get_measure(source, measure, material, medium):

    scatterer = SphereSet(
        diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
        material=material,
        medium=medium
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
    experiment.get(measure, as_numpy=True)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])