"""
Parameter-grid scaling
=======================

Measure the effect of expanding two sweep dimensions: particle diameter and
material refractive index.
"""

import time

import matplotlib.pyplot as plt
import numpy as np

from PyMieSim import Experiment, GaussianSet, PolarizationSet, SphereSet, ureg


def run_grid(size: int) -> tuple[float, int]:
    source = GaussianSet(
        wavelength=[600] * ureg.nanometer,
        polarization=PolarizationSet(angles=[0] * ureg.degree),
        optical_power=[1e-3] * ureg.watt,
        numerical_aperture=[0.2],
    )
    scatterer = SphereSet(
        diameter=np.linspace(100, 1000, size) * ureg.nanometer,
        material=np.linspace(1.3, 1.8, size),
        medium=[1.0],
    )
    experiment = Experiment(scatterer_set=scatterer, source_set=source)

    start = time.perf_counter()
    result = experiment.get("Qsca", as_numpy=True)
    return time.perf_counter() - start, int(np.asarray(result).size)


grid_sizes = [2, 4, 8, 12]
measurements = [run_grid(size) for size in grid_sizes]
runtime_seconds, result_sizes = np.asarray(measurements).T
configuration_counts = np.square(grid_sizes)

figure, axis = plt.subplots()
axis.plot(configuration_counts, runtime_seconds, marker="o")
axis.set(
    xlabel="Number of parameter combinations",
    ylabel="Runtime [s]",
    title="PyMieSim parameter-grid scaling",
)
axis.grid(True, alpha=0.3)
figure.tight_layout()

for combinations, runtime, results in zip(
    configuration_counts, runtime_seconds, result_sizes
):
    print(f"{combinations:4d} combinations | {runtime:.4f} s | {int(results)} results")

plt.show()
