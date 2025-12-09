import pytest
import numpy
from TypedUnit import ureg

from PyMieSim.experiment.scatterer import Sphere, CoreShell
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment.setup import Setup


@pytest.fixture
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,  # 750 nm
        polarization=30 * ureg.degree,  # Polarization in ureg.degrees
        optical_power=1 * ureg.watt,  # Power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical Aperture
    )


@pytest.mark.parametrize("metric", ["Qsca", "Cext", "Qabs", "Qpr", "g"])
def test_no_shell(metric, source):
    diameter = numpy.linspace(100, 1_000, 20) * ureg.nanometer

    sphere = Sphere(
        diameter=diameter,
        property=1.5 * ureg.RIU,
        source=source,
        medium_property=1.33 * ureg.RIU,
    )

    coreshell = CoreShell(
        core_diameter=diameter,
        shell_thickness=0 * ureg.nanometer,
        core_property=1.5 * ureg.RIU,
        shell_property=1.6 * ureg.RIU,
        source=source,
        medium_property=1.33 * ureg.RIU,
    )

    experiment_sphere = Setup(scatterer=sphere, source=source)
    experiment_core_shell = Setup(scatterer=coreshell, source=source)

    data_sphere = experiment_sphere.get(metric, add_units=False).to_numpy().squeeze()
    data_coreshell = (
        experiment_core_shell.get(metric, add_units=False).to_numpy().squeeze()
    )

    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a Sphere and a (no-shell) CoreShell."


@pytest.mark.parametrize("metric", ["Qsca", "Qext", "Qabs", "Qpr", "g"])
def test_shell_equal_core(metric, source):
    diameter_sphere = numpy.linspace(200, 1_000, 20) * ureg.nanometer
    diameter_core = numpy.linspace(100, 900, 20) * ureg.nanometer

    sphere = Sphere(
        diameter=diameter_sphere,
        property=1.5 * ureg.RIU,
        source=source,
        medium_property=1.33 * ureg.RIU,
    )

    coreshell = CoreShell(
        core_diameter=diameter_core,
        shell_thickness=100 * ureg.nanometer,
        core_property=1.5 * ureg.RIU,
        shell_property=1.5 * ureg.RIU,
        source=source,
        medium_property=1.33 * ureg.RIU,
    )

    experiment_sphere = Setup(scatterer=sphere, source=source)
    experiment_core_shell = Setup(scatterer=coreshell, source=source)

    data_sphere = experiment_sphere.get(metric, add_units=False).to_numpy().squeeze()
    data_coreshell = (
        experiment_core_shell.get(metric, add_units=False).to_numpy().squeeze()
    )

    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a Sphere and a (no-shell) CoreShell."


@pytest.mark.parametrize("metric", ["Qsca", "Qext", "Qabs", "Qpr", "g"])
def test_only_shell(metric, source):
    diameter_sphere = numpy.linspace(200, 1_000, 20) * ureg.nanometer
    diameter_shell = numpy.linspace(100, 900, 20) * ureg.nanometer

    sphere = Sphere(
        diameter=diameter_sphere,
        property=1.5 * ureg.RIU,
        source=source,
        medium_property=1.33 * ureg.RIU,
    )

    coreshell = CoreShell(
        core_diameter=100 * ureg.nanometer,
        shell_thickness=diameter_shell,
        core_property=1.5 * ureg.RIU,
        shell_property=1.5 * ureg.RIU,
        source=source,
        medium_property=1.33 * ureg.RIU,
    )

    experiment_sphere = Setup(scatterer=sphere, source=source)
    experiment_core_shell = Setup(scatterer=coreshell, source=source)

    data_sphere = experiment_sphere.get(metric, add_units=False).to_numpy().squeeze()
    data_coreshell = (
        experiment_core_shell.get(metric, add_units=False).to_numpy().squeeze()
    )

    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a Sphere and a (no-shell) CoreShell."


@pytest.mark.parametrize("metric", ["Csca", "Cext", "Cabs"])
def test_shell_is_medium(metric, source):
    diameter = numpy.linspace(400, 800, 100) * ureg.nanometer

    sphere = Sphere(
        diameter=diameter,
        property=1.5 * ureg.RIU,
        source=source,
        medium_property=1.401 * ureg.RIU,
    )

    coreshell = CoreShell(
        core_diameter=diameter,
        shell_thickness=500 * ureg.nanometer,
        core_property=1.5 * ureg.RIU,
        shell_property=1.4 * ureg.RIU,
        source=source,
        medium_property=1.401 * ureg.RIU,
    )

    experiment_sphere = Setup(scatterer=sphere, source=source)
    experiment_core_shell = Setup(scatterer=coreshell, source=source)

    data_sphere = experiment_sphere.get(metric, add_units=False).to_numpy().squeeze()
    data_coreshell = (
        experiment_core_shell.get(metric, add_units=False).to_numpy().squeeze()
    )

    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a Sphere and a (no-shell) CoreShell."


if __name__ == "__main__":
    pytest.main(["-W error", "-s", __file__])
