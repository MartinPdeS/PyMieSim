"""Tests for the stable Python-facing API."""

from PyMieSim import Measure, PlaneWave, PolarizationState, Simulation, Sphere, SimulationResult, SimulationResults, Qsca, ureg


def test_public_api_builds_a_simulation():
    simulation = Simulation(
        scatterer=Sphere(
            diameter=100 * ureg.nanometer,
            material=1.5,
            medium=1.0,
        ),
        source=PlaneWave(
            amplitude=1 * ureg.volt / ureg.meter,
            polarization=PolarizationState(angle=0 * ureg.degree),
            wavelength=1550 * ureg.nanometer,
        ),
    )

    assert Qsca == Measure.QSCA.value
    assert "Qsca" in simulation.available_measures
    assert "Simulation" in repr(simulation)
    assert simulation.run("Qsca") is not None


def test_typed_single_result_preserves_units():
    simulation = Simulation(
        scatterer=Sphere(diameter=100 * ureg.nanometer, material=1.5, medium=1.0),
        source=PlaneWave(
            amplitude=1 * ureg.volt / ureg.meter,
            polarization=PolarizationState(angle=0 * ureg.degree),
            wavelength=1550 * ureg.nanometer,
        ),
    )

    result = simulation.run(Measure.QSCA, as_result=True)

    assert isinstance(result, SimulationResult)
    assert result.measure == "Qsca"
    assert result.units == ureg.dimensionless


def test_typed_single_results_are_mapping_for_multiple_measures():
    simulation = Simulation(
        scatterer=Sphere(diameter=100 * ureg.nanometer, material=1.5, medium=1.0),
        source=PlaneWave(
            amplitude=1 * ureg.volt / ureg.meter,
            polarization=PolarizationState(angle=0 * ureg.degree),
            wavelength=1550 * ureg.nanometer,
        ),
    )

    results = simulation.run(Measure.QSCA, Measure.QEXT, as_result=True)

    assert isinstance(results, SimulationResults)
    assert tuple(results) == ("Qsca", "Qext")
