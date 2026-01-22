#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyOptik import Material
from PyMieSim.units import ureg

from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import InfiniteCylinder
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment import Setup

# Configure the medium and core materials for the sphere
properties = [Material.silver, Material.fused_silica, 1.4 * ureg.RIU]
medium_properties = [Material.water, 1.1 * ureg.RIU]

# Measures to be tested
measures = InfiniteCylinder.available_measure_list

gaussian_source = Gaussian(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)

planewave_source = PlaneWave(
    wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
    polarization=0 * ureg.degree,
    amplitude=1 * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize(
    "medium_refractive_index", medium_properties, ids=[f"Medium:{m}" for m in medium_properties]
)
@pytest.mark.parametrize(
    "source", sources, ids=[f"Source:{m.__class__.__name__}" for m in sources]
)
@pytest.mark.parametrize(
    "refractive_index", properties, ids=[f"Property:{m}" for m in properties]
)
@pytest.mark.parametrize("measure", measures)
def test_measure(measure, source, medium_refractive_index, refractive_index):
    # Configure the cylindrical scatterer
    scatterer = InfiniteCylinder(
        diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
        source=source,
        medium_refractive_index=medium_refractive_index,
        refractive_index=refractive_index,
    )

    # Configure the detector
    detector = Photodiode(
        NA=0.2 * ureg.AU,
        polarization_filter=None,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100 * ureg.AU,
    )

    # Configure and run the experiment
    experiment = Setup(scatterer=scatterer, source=source, detector=detector)

    # Execute measurement
    experiment.get(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
