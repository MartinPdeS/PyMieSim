#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave
from PyMieSim.detector import Photodiode

samplings = [100, 200, 300, 400]


@pytest.mark.parametrize('sampling', samplings)
def test_photodiode(sampling: int):
    source = PlaneWave(
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

    detector = Photodiode(
        NA=0.2,
        sampling=sampling,
        gamma_offset=0,
        phi_offset=0
    )

    detector.get_footprint(scatterer=scatterer)

# -
