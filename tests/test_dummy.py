#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

import numpy


def test_simple():
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
        index=1.5,
        medium_index=1.0
    )

    phi, theta = numpy.mgrid[
        -0.3: 0.3: 100, 0: numpy.pi: 100
    ]

    phi = phi.ravel() + numpy.pi / 2

    theta = theta.ravel()

    a = scatterer.binding.get_fields(phi=phi, theta=theta, r=1.0)

    print(a)


if __name__ == "__main__":
    pytest.main()

# -
