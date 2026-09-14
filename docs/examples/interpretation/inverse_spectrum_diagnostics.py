"""
Diagnose a sphere fit across several wavelengths
================================================

A single scattering measurement can leave diameter and refractive index
ambiguous. A spectrum supplies more structure, while residuals show whether
the fitted isolated-sphere model explains every channel consistently.

The synthetic observation below contains one deliberately biased wavelength.
That channel is assigned a larger uncertainty, so it contributes less to the
weighted objective. ``FitResult.residuals`` uses
``(prediction - observation) / uncertainty``; a point above the model therefore
has a negative standardized residual.
"""
import matplotlib.pyplot as plt
import numpy as np

from PyMieSim import Observation, Parameter, PlaneWave, PolarizationState, Simulation, Sphere, fit_parameters, ureg

wavelengths = np.linspace(450, 750, 9) * ureg.nanometer


def scattering_spectrum(diameter, refractive_index):
    """Return isolated-sphere scattering efficiencies at all wavelengths."""
    spectrum = np.empty(wavelengths.size)
    for index, wavelength in enumerate(wavelengths):
        source = PlaneWave(
            wavelength=wavelength,
            polarization=PolarizationState(angle=0 * ureg.degree),
            amplitude=1 * ureg.volt / ureg.meter,
        )
        simulation = Simulation(
            scatterer=Sphere(
                diameter=diameter,
                material=refractive_index,
                medium=1.0,
            ),
            source=source,
        )
        spectrum[index] = simulation.get("Qsca").magnitude
    return spectrum


observed = scattering_spectrum(
    diameter=310 * ureg.nanometer,
    refractive_index=1.58,
)
uncertainty = np.full_like(observed, 0.015)
observed[5] += 0.08
uncertainty[5] = 0.08


def model(parameters):
    return scattering_spectrum(
        diameter=parameters["diameter"],
        refractive_index=parameters["refractive_index"],
    )


result = fit_parameters(
    model=model,
    observation=Observation(values=observed, uncertainty=uncertainty),
    parameters=[
        Parameter(
            name="diameter",
            initial=240 * ureg.nanometer,
            bounds=(150 * ureg.nanometer, 450 * ureg.nanometer),
        ),
        Parameter(name="refractive_index", initial=1.45, bounds=(1.2, 2.0)),
    ],
)

wavelength_values = wavelengths.to("nanometer").magnitude
figure, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True, height_ratios=(2, 1))
axes[0].errorbar(wavelength_values, observed, yerr=uncertainty, fmt="o", label="Synthetic observation")
axes[0].plot(wavelength_values, result.prediction, linewidth=2, label="Fitted isolated-sphere model")
axes[0].set(ylabel="Qsca", title="Multi-wavelength fitting constrains a spectral response")
axes[0].legend()
axes[1].axhline(0, color="black", linewidth=1)
axes[1].plot(wavelength_values, result.residuals, "o-")
axes[1].set(xlabel="Vacuum wavelength [nm]", ylabel="Standardized residual")
for axis in axes:
    axis.grid(alpha=0.2)
plt.show()
