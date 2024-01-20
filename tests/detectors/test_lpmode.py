#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import LPMode

mode_numbers = [
    "LP01",
    "LP11",
    "LP21"
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
        n_medium=1.0
    )

    detector = LPMode(
        mode_number=mode_number,
        NA=0.2,
        sampling=100,
        gamma_offset=0,
        phi_offset=0
    )

    detector.get_footprint(scatterer=scatterer)

# -
