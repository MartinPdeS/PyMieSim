"""
Runtime scaling with sweep size
================================

Measure how long a ``Qsca`` calculation takes as the number of sphere
diameters increases.
"""

import time

import matplotlib.pyplot as plt
import numpy as np

from PyMieSim import Experiment, GaussianSet, PolarizationSet, SphereSet, ureg


def run_sweep(number_of_diameters: int) -> tuple[float, int]:
    source = GaussianSet(
        wavelength=[600] * ureg.nanometer,
        polarization=PolarizationSet(angles=[0] * ureg.degree),
        optical_power=[1e-3] * ureg.watt,
        numerical_aperture=[0.2],
    )
    scatterer = SphereSet(
        diameter=np.linspace(100, 1000, number_of_diameters) * ureg.nanometer,
        material=[1.5],
        medium=[1.0],
    )
    experiment = Experiment(scatterer_set=scatterer, source_set=source)

    start = time.perf_counter()
    values = experiment.get("Qsca").as_numpy()
    elapsed = time.perf_counter() - start
    return elapsed, int(np.asarray(values).size)


sweep_sizes = [10, 25, 50, 100, 200, 400]
measurements = [run_sweep(size) for size in sweep_sizes]
runtime_seconds, result_sizes = np.asarray(measurements).T

for size, runtime, result_size in zip(sweep_sizes, runtime_seconds, result_sizes):
    print(f"{size:4d} diameters | {runtime:.4f} s | {int(result_size)} results")

figure, axis = plt.subplots()
axis.plot(sweep_sizes, runtime_seconds, marker="o")
axis.set(
    xlabel="Number of diameters",
    ylabel="Runtime [s]",
    title="PyMieSim runtime scaling",
)
axis.grid(True, alpha=0.3)
figure.tight_layout()
plt.show()
