#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian


scott_prahl_values = {
    'Qsca': 1.1759,
    'Qext': 2.6257,
    'g': 0.80335,
}


@pytest.mark.parametrize('measure_str', scott_prahl_values.keys(), ids=scott_prahl_values.keys())
def test_validation_scott_prahl(measure_str):
    source = Gaussian(
        wavelength=1e-6,
        polarization=0,
        optical_power=1,
        NA=0.3
    )

    scatterer = Sphere(
        diameter=1e-6,
        index=1.5 + 0.5j,
        source=source,
        medium_index=1.0
    )

    scott_prahl_value = scott_prahl_values.get(measure_str)
    pymiesim_value = getattr(scatterer, measure_str)

    if not numpy.isclose(pymiesim_value, scott_prahl_value, atol=0, rtol=1e-3):
        raise ValueError('Mismatch with testing values')


if __name__ == "__main__":
    pytest.main()

# -
