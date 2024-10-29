import pytest
import numpy
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU


size = 10
ones = numpy.ones(size)
wavelength = numpy.linspace(100, 400, size) * nanometer
medium_property = 1.0 * RIU
property = 1.5 * RIU
diameter = 100 * nanometer
polarization = 0 * degree
optical_power = 200 * watt
NA = 0.2 * AU


@pytest.fixture
def sequential_experiment():

    source = Gaussian(
        wavelength=wavelength,
        polarization=ones * polarization,
        optical_power=ones * optical_power,
        NA=ones * NA
    )

    scatterer = Sphere(
        diameter=ones * diameter,
        property=ones * property,
        medium_property=ones * medium_property,
        source=source
    )

    return Setup(scatterer=scatterer, source=source)


@pytest.fixture
def standard_experiment():
    source = Gaussian(
        wavelength=wavelength,
        polarization=polarization,
        optical_power=optical_power,
        NA=NA
    )

    scatterer = Sphere(
        diameter=diameter,
        property=property,
        medium_property=medium_property,
        source=source
    )

    return Setup(scatterer=scatterer, source=source)


def test_sequential_vs_standard(sequential_experiment, standard_experiment):
    values_seq = sequential_experiment.get_sequential('Qsca')
    values_std = standard_experiment.get('Qsca', add_units=False).values.squeeze()

    assert numpy.allclose(values_seq, values_std), "Mismatch between sequential and structured computation"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
