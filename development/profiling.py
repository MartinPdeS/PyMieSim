import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU
import cProfile

from pathlib import Path

source = Gaussian(
    wavelength=400 * nanometer,
    polarization=[0] * degree,
    optical_power=1e-6 * watt,
    NA=0.2 * AU
)

scatterer = Sphere(
    diameter=np.linspace(300, 1000, 1000) * nanometer,
    property=[1.2, 1.25] * RIU,
    medium_property=[1.0] * RIU,
    source=source
)

experiment = Setup(scatterer=scatterer, source=source)

profiler = cProfile.Profile()


profiler.enable()

dataframe = experiment.get('a1')

# Stop profiling
profiler.disable()

# Dump the profile data to a file
profiler.dump_stats(Path(__file__).parent / "output_file.prof")
