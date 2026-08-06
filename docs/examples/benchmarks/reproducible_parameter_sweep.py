"""
Reproducible parameter-sweep benchmark
======================================

Benchmark a fixed PyMieSim parameter grid while recording environment metadata.
The script is deterministic; runtime values are machine-dependent and should
be compared only under the same environment.
"""

from __future__ import annotations

import platform
import statistics
import time

import numpy as np

from PyMieSim import Experiment, GaussianSet, PolarizationSet, SphereSet, ureg


def build_experiment(size: int = 100) -> Experiment:
    source = GaussianSet(
        wavelength=[600] * ureg.nanometer,
        polarization=PolarizationSet(angles=[0] * ureg.degree),
        optical_power=[1e-3] * ureg.watt,
        numerical_aperture=[0.2],
    )
    scatterer = SphereSet(
        diameter=np.linspace(100, 1000, size) * ureg.nanometer,
        material=[1.5],
        medium=[1.0],
    )
    return Experiment(scatterer_set=scatterer, source_set=source)


experiment = build_experiment()
for _ in range(2):
    experiment.get("Qsca").as_numpy()

timings = []
for _ in range(5):
    start = time.perf_counter()
    values = experiment.get("Qsca").as_numpy()
    timings.append(time.perf_counter() - start)

print(f"Python: {platform.python_version()}")
print(f"Platform: {platform.platform()}")
print(f"PyMieSim grid: {experiment.array_shape}")
print(f"Results: {np.asarray(values).size}")
print(f"Median runtime: {statistics.median(timings):.6f} s")
print(f"All timings: {[round(value, 6) for value in timings]}")
