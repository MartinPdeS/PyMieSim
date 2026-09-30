"""
Fit a scalar scale factor
=========================

Fit one bounded parameter to a small analytic forward model. This is the
smallest complete inverse-problem setup: a callable, an observation, and one
or more ``Parameter`` objects.
"""
import numpy as np

from PyMieSim import Observation, Parameter, fit_parameters

x = np.array([1.0, 2.0, 3.0])
observation = Observation(values=np.array([2.0, 4.0, 6.0]))


def model(parameters):
    return parameters["scale"] * x


result = fit_parameters(
    model=model,
    observation=observation,
    parameters=[Parameter(name="scale", initial=0.8, bounds=(0.0, 4.0))],
)

print(result.summary())
