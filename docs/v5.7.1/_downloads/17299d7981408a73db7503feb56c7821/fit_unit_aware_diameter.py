"""
Fit a unit-aware diameter
=========================

Parameters can carry Pint units. The model receives values with the original
units, so the physics code does not need to know how the optimizer stores its
internal scalar coordinates.
"""
import numpy as np

from PyMieSim import Observation, Parameter, fit_parameters, ureg

observation = Observation(values=240 * ureg.nanometer)


def model(parameters):
    diameter = parameters["diameter"]
    return 2 * diameter


result = fit_parameters(
    model=model,
    observation=observation,
    parameters=[
        Parameter(
            name="diameter",
            initial=150 * ureg.nanometer,
            bounds=(50 * ureg.nanometer, 300 * ureg.nanometer),
        )
    ],
)

print(result.parameters["diameter"])
