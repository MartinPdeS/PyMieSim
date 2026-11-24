#!/usr/bin/env python
# -*- coding: utf-8 -*-

import PyOptik
from PyMieSim import from_pyoptik

import pytest
import numpy as np

from PyMieSim.units import ureg
from PyMieSim.experiment.detector_set import PhotodiodeSet
from PyMieSim.experiment.scatterer_set import InfiniteCylinderSet
from PyMieSim.experiment.source_set import GaussianSet, PlaneWaveSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim import Experiment


properties = [
    [from_pyoptik(PyOptik.load_material(name="main/Ag/Johnson"))],
    [from_pyoptik(PyOptik.load_material(name="main/SiO2/Malitson"))],
    [1.4]
]
medium_properties = [
    [from_pyoptik(PyOptik.load_material(name="main/H2O/Daimon-19.0C"), medium=True)],
    [1.1]
]

measures = InfiniteCylinderSet.available_measure_list


polarization = PolarizationSet(angles=0 * ureg.degree)

gaussian_source = GaussianSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.2],
)

planewave_source = PlaneWaveSet(
    wavelength=np.linspace(600, 1000, 15) * ureg.nanometer,
    polarization=polarization,
    amplitude=[1] * ureg.volt / ureg.meter,
)

sources = [gaussian_source, planewave_source]


@pytest.mark.parametrize( "medium", medium_properties,
    ids=[f"Medium:{m}" for m in medium_properties],
)
@pytest.mark.parametrize(
    "source",
    sources,
    ids=[f"Source:{m.__class__.__name__}" for m in sources],
)
@pytest.mark.parametrize( "material", properties,
    ids=[f"Property:{m}" for m in properties],
)
@pytest.mark.parametrize("measure", measures)
def test_measure(measure, source, medium, material):
    scatterer = InfiniteCylinderSet(
        diameter=np.linspace(400, 1400, 10) * ureg.nanometer,
        material=material,
        medium=medium
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
        detector_set=detector,
    )

    values = experiment.run(measure, drop_unique_level=True).as_numpy()
    assert values.shape == (15, 10)
    assert np.isfinite(values).all()


@pytest.mark.parametrize("measure", ["Qsca", "Qext"])
def test_cylinder_efficiencies_are_invariant_at_fixed_size_parameter_and_index_ratio(measure):
    """Match both incident polarizations in a dielectric host and in vacuum."""
    vacuum_wavelength = 600 * ureg.nanometer
    diameters = [180, 400] * ureg.nanometer
    material_index = 1.6 + 0.05j
    host_index = 1.33
    polarization = PolarizationSet(angles=[0, 90] * ureg.degree)

    embedded = Experiment(
        scatterer_set=InfiniteCylinderSet(diameter=diameters, material=material_index, medium=host_index),
        source_set=PlaneWaveSet(
            wavelength=vacuum_wavelength, polarization=polarization, amplitude=1 * ureg.volt / ureg.meter
        ),
    )
    equivalent = Experiment(
        scatterer_set=InfiniteCylinderSet(diameter=diameters, material=material_index / host_index, medium=1.0),
        source_set=PlaneWaveSet(
            wavelength=vacuum_wavelength / host_index, polarization=polarization, amplitude=1 * ureg.volt / ureg.meter
        ),
    )

    np.testing.assert_allclose(embedded.run(measure).as_numpy(), equivalent.run(measure).as_numpy(), rtol=1e-12, atol=0)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])