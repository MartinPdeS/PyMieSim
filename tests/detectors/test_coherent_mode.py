#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import CoherentMode

mode_numbers = [
    "LP01:00",
    "LP11:45",
    "LP21:90",
    "LG01:00",
    "LG11:45",
    "LG21:90",
    "HG01:00",
    "HG11:45",
    "HG21:90"
]


@pytest.mark.parametrize('mode_number', mode_numbers)
def test_lp_modes(mode_number: str):
    source = Gaussian(
        wavelength=750e-9,
        polarization_value=0,
        polarization_type='linear',
        optical_power=1,
        NA=0.3
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        index=1.4,
        medium_index=1.0
    )

    detector = CoherentMode(
        mode_number=mode_number,
        NA=0.2,
        sampling=100,
        gamma_offset=0,
        phi_offset=0
    )

    detector.get_footprint(scatterer=scatterer)

    detector.plot()

# -
