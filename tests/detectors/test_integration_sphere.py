#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import IntegratingSphere

samplings = [100]


@pytest.mark.parametrize('sampling', samplings)
def test_photodiode(sampling: int):
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

    detector = IntegratingSphere(sampling=sampling)

    detector.get_footprint(scatterer=scatterer)

    detector.plot()

# -
