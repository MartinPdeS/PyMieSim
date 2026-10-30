#!/usr/bin/env python
# -*- coding: utf-8 -*-

import PyOptik
from PyMieSim import from_pyoptik

import pytest
import numpy as np

from PyMieSim.units import ureg
from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.scatterer_set import CoreShellSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim import Experiment

core_properties = [
    [from_pyoptik(PyOptik.load_material(name="main/Ag/Johnson"))],
    [from_pyoptik(PyOptik.load_material(name="main/SiO2/Malitson"))],
    [1.4]
]
shell_properties = [
    [from_pyoptik(PyOptik.load_material(name="main/Ag/Johnson"))],
    [from_pyoptik(PyOptik.load_material(name="main/SiO2/Malitson"))],
    [1.4]
]
medium_properties = [
    [from_pyoptik(PyOptik.load_material(name="main/H2O/Daimon-19.0C"), medium=True)],
    [1.1]
]


polarization_set = PolarizationSet(
    angles=[0] * ureg.radian
)


gaussian_source = GaussianSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

planewave_source = PlaneWaveSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization_set,
    amplitude=[1] * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize(
    "medium",
    medium_properties,
    ids=[f"Medium:{m}" for m in medium_properties]
)

@pytest.mark.parametrize(
    "core_material",
    core_properties,
    ids=[f"Core:{m}" for m in core_properties]
)

@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources]
)

@pytest.mark.parametrize(
    "shell_material",
    shell_properties,
    ids=[f"Shell:{m}" for m in shell_properties]
)

@pytest.mark.parametrize(
    "measure",
    CoreShellSet.available_measure_list[1:]
)

def test_measure(measure, source, core_material, shell_material, medium):
    kwargs = {}

    scatterer = CoreShellSet(
        core_diameter=np.linspace(800, 1000, 10) * ureg.nanometer,
        shell_thickness=[300] * ureg.nanometer,
        core_material=core_material,
        shell_material=shell_material,
        medium=medium,
        **kwargs
    )

    detector = PhotodiodeSet(
        numerical_aperture=[0.2],
        gamma_offset=[0] * ureg.degree,
        phi_offset=[0] * ureg.degree,
        sampling=[100],
    )

    experiment = Experiment(
        scatterer_set=scatterer,
        source_set=source,
        detector_set=detector
    )

    experiment.run(measure, drop_unique_level=False)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])