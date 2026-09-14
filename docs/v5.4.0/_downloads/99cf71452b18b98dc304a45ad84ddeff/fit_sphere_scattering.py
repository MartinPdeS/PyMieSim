"""
Fit an isolated-sphere scattering measurement
==============================================

Recover a sphere diameter from one dimensionless scattering-efficiency
measurement. This uses the isolated-particle Mie model. It does not represent
multiple scattering, near-field coupling, or positional correlations in dense
particle systems.
"""
from PyMieSim import (
    Observation, Parameter, PlaneWave, PolarizationState, Simulation, Sphere,
    fit_parameters, ureg,
)

source = PlaneWave(
    wavelength=532 * ureg.nanometer,
    polarization=PolarizationState(angle=0 * ureg.degree),
    amplitude=1 * ureg.volt / ureg.meter,
)
target_diameter = 180 * ureg.nanometer
target = Simulation(
    scatterer=Sphere(diameter=target_diameter, material=1.59, medium=1.0),
    source=source,
).get("Qsca")


def model(parameters):
    sphere = Sphere(diameter=parameters["diameter"], material=1.59, medium=1.0)
    return Simulation(scatterer=sphere, source=source).get("Qsca")


result = fit_parameters(
    model=model,
    observation=Observation(values=target),
    parameters=[
        Parameter(
            name="diameter",
            initial=130 * ureg.nanometer,
            bounds=(80 * ureg.nanometer, 260 * ureg.nanometer),
        )
    ],
    max_iterations=250,
    show_progress=True,
)

print(result.summary())
