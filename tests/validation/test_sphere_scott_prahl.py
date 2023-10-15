#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy

from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave


scott_prahl_values = {
    'Qsca': 1.1759,
    'Qext': 2.6257,
    'g': 0.80335,
}


@pytest.mark.parametrize('measure_str', scott_prahl_values.keys(), ids=scott_prahl_values.keys())
def test_validation_scott_prahl(measure_str):
    source = PlaneWave(
        wavelength=1e-6,
        polarization=0,
        amplitude=1
    )

    scatterer = Sphere(
        diameter=1e-6,
        index=1.5 + 0.5j,
        source=source
    )

    scott_prahl_value = scott_prahl_values.get(measure_str)
    pymiesim_value = getattr(scatterer, measure_str)

    if not numpy.isclose(pymiesim_value, scott_prahl_value, atol=0, rtol=1e-3):
        raise ValueError('Mismatch with testing values')


# -
