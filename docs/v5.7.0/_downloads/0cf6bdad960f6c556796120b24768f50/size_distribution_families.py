"""
Compare size distribution families
=================================

Number-based size populations averaged in the non-interacting approximation.
No interparticle coupling, multiple scattering, or interparticle interference
is included. Refer to PackLab for correlation-based dependent-scattering
approximations when particle correlations matter.
"""
import matplotlib.pyplot as plt
import numpy as np

from PyMieSim import Experiment, ParticleSizeDistribution, PlaneWaveSet, PolarizationSet, SphereSet, ureg

sizes = ParticleSizeDistribution
populations = {
    "Monodisperse": sizes.monodisperse(200 * ureg.nanometer),
    "Uniform": sizes.uniform(100 * ureg.nanometer, 300 * ureg.nanometer),
    "Lognormal": sizes.lognormal(200 * ureg.nanometer, 1.2),
    "Truncated normal": sizes.truncated_normal(
        200 * ureg.nanometer, 40 * ureg.nanometer,
        minimum_diameter=100 * ureg.nanometer, maximum_diameter=300 * ureg.nanometer,
    ),
    "Triangular": sizes.triangular(100 * ureg.nanometer, 200 * ureg.nanometer, 300 * ureg.nanometer),
}
populations["Mixture (75% / 25%)"] = sizes.mixture(
    [sizes.monodisperse(100 * ureg.nanometer), populations["Uniform"]], [3, 1],
)
source = PlaneWaveSet(
    wavelength=np.linspace(400, 800, 40) * ureg.nanometer,
    polarization=PolarizationSet(angles=0 * ureg.degree),
    amplitude=[1] * ureg.volt / ureg.meter,
)
figure, axis = plt.subplots()
for label, distribution in populations.items():
    experiment = Experiment(
        source_set=source,
        scatterer_set=SphereSet(diameter=distribution.diameters, material=[1.5], medium=[1.0]),
    )
    result = experiment.average_size_distribution(distribution, "Csca", as_result=True).to("nanometer ** 2")
    axis.plot(result.coords["source:wavelength"] * 1e9, result.magnitude, label=label)
axis.set(xlabel="Wavelength [nm]", ylabel="Mean scattering cross section per particle [nm²]",
         title="Size distributions — non-interacting particles")
axis.legend()
figure.tight_layout()
plt.show()
