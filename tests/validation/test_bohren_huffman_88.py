#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from PyMieSim.units import ureg


# PyMieSim imports
from PyMieSim.directories import validation_data_path
from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup


def test_validation():

    theoretical_data = np.genfromtxt(
        f"{validation_data_path}/bohren_huffman/figure_88.csv", delimiter=","
    ) / (ureg.meter * ureg.meter)

    wavelength = 632.8 * ureg.nanometer
    polarization_values = PolarizationSet(angles=[0, 90] * ureg.degree)
    optical_power = 1e-3 * ureg.watt
    NA = 0.2
    diameters = np.geomspace(10, 6000, 800) * ureg.nanometer
    index = 1.55
    medium_index = 1.335

    volumes = np.pi * (diameters / 2) ** 2

    source = GaussianSet(
        wavelength=wavelength,
        polarization=polarization_values,
        optical_power=optical_power,
        numerical_aperture=NA,
    )

    scatterer = InfiniteCylinderSet(
        diameter=diameters, material=index, medium=medium_index
    )

    experiment = Setup(scatterer_set=scatterer, source_set=source)

    csca_data =  experiment.get("Csca", as_numpy=True)

    normalized_csca = csca_data / volumes.to_base_units() * 1e-4 / 100

    assert np.allclose(normalized_csca[0], theoretical_data[0], atol=0, rtol=1e-3), (
        f"Mismatch in Csca for polarization 0: PyMieSim value = {normalized_csca[0]}, "
        f"Theoretical value = {theoretical_data[0]}"
    )
    assert np.allclose(normalized_csca[1], theoretical_data[1], atol=0, rtol=1e-3), (
        f"Mismatch in Csca for polarization 90: PyMieSim value = {normalized_csca[1]}, "
        f"Theoretical value = {theoretical_data[1]}"
    )


if __name__ == "__main__":
    pytest.main(["-W error", "-s", __file__])
