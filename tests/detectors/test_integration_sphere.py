#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave
from PyMieSim.detector import IntegratingSphere

samplings = [100]


@pytest.mark.parametrize('sampling', samplings)
def test_photodiode(sampling: int):
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

    detector = IntegratingSphere(sampling=sampling)

    detector.get_footprint(scatterer=scatterer)

    detector.plot()

# -
