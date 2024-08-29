#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from unittest.mock import patch

@patch('matplotlib.pyplot.show')
def test_sphere_plottings(mock_show_plt):
    source = Gaussian(
        wavelength=750e-9,
        polarization=0,
        optical_power=1,
        NA=0.3
    )
    scatterer = Sphere(
        diameter=100e-9,
        source=source,
        index=1.4,
        medium_index=1.1
    )

    detector = Photodiode(
        NA=0.1,
        phi_offset=0,
        gamma_offset=0,
        polarization_filter=0
    )

    data = scatterer.get_footprint(detector=detector)
    data.plot()
    assert data is not None, "Plotting data should not be None"


if __name__ == "__main__":
    pytest.main([__file__])


# -
