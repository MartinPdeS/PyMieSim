"""
Energy-conservation validation
==============================

Check the identity ``Qext = Qsca + Qabs`` over a wavelength sweep for a
lossy sphere. This is an internal consistency validation that is independent
of plotting and can be reused as a compact regression example.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyMieSim import Gaussian, Measure, PolarizationState, Simulation, Sphere, ureg


wavelengths = np.linspace(450, 750, 25) * ureg.nanometer
residuals = []

for wavelength in wavelengths:
    source = Gaussian(
        wavelength=wavelength,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1e-3 * ureg.watt,
        numerical_aperture=0.2,
    )
    simulation = Simulation(
        scatterer=Sphere(
            diameter=180 * ureg.nanometer,
            material=1.6 + 0.03j,
            medium=1.0,
        ),
        source=source,
    )
    qext = simulation.run(Measure.QEXT).magnitude
    qsca = simulation.run(Measure.QSCA).magnitude
    qabs = simulation.run(Measure.QABS).magnitude
    residuals.append(float(qext - qsca - qabs))

residuals = np.asarray(residuals)
print(f"Maximum absolute residual: {np.max(np.abs(residuals)):.3e}")
assert np.max(np.abs(residuals)) < 1e-10

figure, axis = plt.subplots()
axis.plot(wavelengths.to("nanometer").magnitude, residuals, marker="o")
axis.set(xlabel="Wavelength [nm]", ylabel="Qext - Qsca - Qabs")
axis.axhline(0, color="black", linewidth=0.8)
figure.tight_layout()
plt.close(figure)
