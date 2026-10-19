"""
Fit a scalar scale factor
=========================

Fit one bounded parameter to a small analytic forward model. This is the
smallest complete inverse-problem setup: a callable, an observation, and one
or more ``FitParameter`` objects.

The fitted ``scale`` is approximately 2.0 (dimensionless). Predictions and
residuals have shape ``(3,)``. With these exact synthetic observations, the
residuals approach zero.
"""
import numpy as np

from PyMieSim import FitParameter, MeasuredData, fit_parameters

x = np.array([1.0, 2.0, 3.0])
observation = MeasuredData(values=np.array([2.0, 4.0, 6.0]))


def model(parameters):
    return parameters["scale"] * x


result = fit_parameters(
    model=model,
    observation=observation,
    parameters=[FitParameter(name="scale", initial=0.8, bounds=(0.0, 4.0))],
)

print(result.summary())
print(result.parameters)
