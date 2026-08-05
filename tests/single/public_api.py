"""Tests for the stable Python-facing API."""

from PyMieSim import PlaneWave, PolarizationState, Simulation, Sphere, ureg


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

    assert simulation.setup is not None
    assert "Simulation" in repr(simulation)
    assert simulation.run("Qsca") is not None
