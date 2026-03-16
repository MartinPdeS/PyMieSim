import pytest
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer_set import SphereSet, CoreShellSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.setup import Setup


@pytest.mark.parametrize("metric", ["Qsca", "Cext", "Qabs", "Qpr", "g"])
def test_no_shell(metric):
    diameter = numpy.linspace(100, 1_000, 20) * ureg.nanometer

    sphere = SphereSet(
        diameter=diameter,
        material=1.5,
        medium=1.33,
    )

    coreshell = CoreShellSet(
        core_diameter=diameter,
        shell_thickness=0 * ureg.nanometer,
        core_material=1.5,
        shell_material=1.6,
        medium=1.33,
    )

    source = GaussianSet(
        wavelength=[750] * ureg.nanometer,
        polarization=PolarizationSet(angles=30 * ureg.degree),
        optical_power=[1] * ureg.watt,
        numerical_aperture=[0.3],
    )

    experiment_sphere = Setup(scatterer_set=sphere, source_set=source)
    experiment_core_shell = Setup(scatterer_set=coreshell, source_set=source)

    data_sphere = experiment_sphere.get(metric, as_numpy=True).squeeze()
    data_coreshell = experiment_core_shell.get(metric, as_numpy=True).squeeze()

    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a SphereSet and a (no-shell) CoreShellSet."


@pytest.mark.parametrize("metric", ["Qsca", "Qext", "Qabs", "Qpr", "g"])
def test_shell_equal_core(metric):
    diameter_sphere = numpy.linspace(200, 1_000, 20) * ureg.nanometer
    diameter_core = numpy.linspace(100, 900, 20) * ureg.nanometer

    source = GaussianSet(
        wavelength=[750] * ureg.nanometer,
        polarization=PolarizationSet(angles=30 * ureg.degree),
        optical_power=[1] * ureg.watt,
        numerical_aperture=[0.3],
    )

    sphere = SphereSet(
        diameter=diameter_sphere,
        material=[1.5],
        medium=[1.33],
    )

    coreshell = CoreShellSet(
        core_diameter=diameter_core,
        shell_thickness=[50] * ureg.nanometer,
        core_material=[1.5],
        shell_material=[1.5],
        medium=[1.33],
    )

    experiment_sphere = Setup(scatterer_set=sphere, source_set=source)
    experiment_core_shell = Setup(scatterer_set=coreshell, source_set=source)

    data_sphere = experiment_sphere.get(metric, as_numpy=True).squeeze()
    data_coreshell = experiment_core_shell.get(metric, as_numpy=True).squeeze()


    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a SphereSet and a (no-shell) CoreShellSet."


@pytest.mark.parametrize("metric", ["Qsca", "Qext", "Qabs", "Qpr", "g"])
def test_only_shell(metric):
    diameter_sphere = numpy.linspace(200, 1_000, 20) * ureg.nanometer
    shell_thickness = numpy.linspace(50, 450, 20) * ureg.nanometer

    source = GaussianSet(
        wavelength=[750] * ureg.nanometer,
        polarization=PolarizationSet(angles=30 * ureg.degree),
        optical_power=[1] * ureg.watt,
        numerical_aperture=[0.3],
    )


    sphere = SphereSet(
        diameter=diameter_sphere,
        material=[1.5],
        medium=[1.33],
    )

    coreshell = CoreShellSet(
        core_diameter=[100] * ureg.nanometer,
        shell_thickness=shell_thickness,
        core_material=[1.5],
        shell_material=[1.5],
        medium=[1.33],
    )

    experiment_sphere = Setup(scatterer_set=sphere, source_set=source)
    experiment_core_shell = Setup(scatterer_set=coreshell, source_set=source)

    data_sphere = experiment_sphere.get(metric, as_numpy=True).squeeze()
    data_coreshell = experiment_core_shell.get(metric, as_numpy=True).squeeze()


    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a SphereSet and a (no-shell) CoreShellSet."


@pytest.mark.parametrize("metric", ["Csca", "Cext", "Csca"])
def test_shell_is_medium(metric):
    source = GaussianSet(
        wavelength=[750] * ureg.nanometer,
        polarization=PolarizationSet(angles=30 * ureg.degree),
        optical_power=[1] * ureg.watt,
        numerical_aperture=[0.3],
    )


    diameter = numpy.linspace(400, 800, 100) * ureg.nanometer

    sphere = SphereSet(
        diameter=diameter,
        material=[1.5],
        medium=[1.4],
    )

    coreshell = CoreShellSet(
        core_diameter=diameter,
        shell_thickness=[500] * ureg.nanometer,
        core_material=[1.5],
        shell_material=[1.4],
        medium=[1.4],
    )

    experiment_sphere = Setup(scatterer_set=sphere, source_set=source)
    experiment_core_shell = Setup(scatterer_set=coreshell, source_set=source)

    data_sphere = experiment_sphere.get(metric, as_numpy=True)
    data_coreshell = (
        experiment_core_shell.get(metric, as_numpy=True)
    )

    assert numpy.allclose(
        data_sphere, data_coreshell
    ), "Mismatch between the computed value for a SphereSet and a (no-shell) CoreShellSet."


if __name__ == "__main__":
    pytest.main(["-W error", "-s", __file__])
