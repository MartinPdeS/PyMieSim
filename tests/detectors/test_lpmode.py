#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave
from PyMieSim.detector import LPmode

mode_numbers = [
    "LP01",
    "LP11",
    "LP21"
]


@pytest.mark.parametrize('mode_number', mode_numbers)
def test_lp_modes(mode_number: str):
    source = PlaneWave(
        wavelength=1e-6,
        linear_polarization=0,
        amplitude=1
    )

    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        index=1.4,
        n_medium=1.0
    )

    detector = LPmode(
        mode_number=mode_number,
        NA=0.2,
        sampling=100,
        gamma_offset=0,
        phi_offset=0
    )

    detector.get_footprint(scatterer=scatterer)

# -
