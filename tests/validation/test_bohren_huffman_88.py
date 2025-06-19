#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

# Standard library imports
import numpy as np

# PyMieSim imports
from PyMieSim.directories import validation_data_path
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, watt, AU, RIU, nanometer, meter


def test_validation():

    theoretical_data = np.genfromtxt(f"{validation_data_path}/bohren_huffman/figure_88.csv", delimiter=',') / (meter * meter)

    wavelength = 632.8 * nanometer
    polarization_values = [0, 90] * degree
    optical_power = 1e-3 * watt
    NA = 0.2 * AU
    diameters = np.geomspace(10, 6000, 800) * nanometer
    index = 1.55 * RIU
    medium_index = 1.335 * RIU

    volumes = np.pi * (diameters / 2)**2

    source = Gaussian(
        wavelength=wavelength,
        polarization=polarization_values,
        optical_power=optical_power,
        NA=NA
    )

    scatterer = Cylinder(
        diameter=diameters,
        property=index,
        medium_property=medium_index,
        source=source
    )

    experiment = Setup(scatterer=scatterer, source=source)

    csca_data = experiment.get('Csca', add_units=False).squeeze().values.reshape([-1, diameters.size])
    normalized_csca = csca_data / volumes.to_base_units() * 1e-4 / 100

    assert np.allclose(normalized_csca[0], theoretical_data[0], atol=0, rtol=1e-3), (
        f'Mismatch in Csca for polarization 0: PyMieSim value = {normalized_csca[0]}, '
        f'Theoretical value = {theoretical_data[0]}'
    )
    assert np.allclose(normalized_csca[1], theoretical_data[1], atol=0, rtol=1e-3), (
        f'Mismatch in Csca for polarization 90: PyMieSim value = {normalized_csca[1]}, '
        f'Theoretical value = {theoretical_data[1]}'
    )
if __name__ == "__main__":
    pytest.main(["-W error", "-s", __file__])

