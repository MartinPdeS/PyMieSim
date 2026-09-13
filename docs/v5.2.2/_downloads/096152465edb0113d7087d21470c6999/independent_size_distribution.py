"""
Size distribution: independent scattering only
==============================================

Average isolated-sphere cross sections over a lognormal number distribution.
This is the non-interacting approximation: no interparticle electromagnetic
coupling, multiple scattering, positional correlations, or interference between
particles is modeled. The result is a mean per particle, not a sample total.
"""
import numpy as np
import matplotlib.pyplot as plt

from PyMieSim import (
    Experiment, ParticleSizeDistribution, PlaneWaveSet, PolarizationSet, SphereSet, ureg,
)

source = PlaneWaveSet(
    wavelength=np.linspace(400, 800, 40) * ureg.nanometer,
    polarization=PolarizationSet(angles=0 * ureg.degree),
    amplitude=[1] * ureg.volt / ureg.meter,
)
figure, axis = plt.subplots()
for sampling in (16, 32, 64):
    distribution = ParticleSizeDistribution.lognormal(
        median_diameter=150 * ureg.nanometer,
        geometric_std=1.25,
        sampling=sampling,
    )
    experiment = Experiment(
        source_set=source,
        scatterer_set=SphereSet(diameter=distribution.diameters, material=[1.5], medium=[1.0]),
    )
    result = experiment.average_size_distribution(
        distribution, "Csca", as_result=True,
    ).to("nanometer ** 2")
    axis.plot(result.coords["source:wavelength"] * 1e9, result.magnitude, label=f"{sampling} quadrature nodes"
            )

axis.set(
    xlabel="Wavelength [nm]",
    ylabel="Mean scattering cross section per particle [nm²]",
    title="Lognormal sizes — non-interacting particles only"
)
axis.legend()
figure.tight_layout()
plt.show()
