"""
Sphere: Coupling vs diameter
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from TypedUnit import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material
from PyMieSim.single.mesh import FibonacciMesh

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200 * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=[0.1] * ureg.AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=numpy.linspace(100, 10000, 600) * ureg.nanometer,
    property=Material.BK7,
    medium_property=1.0 * ureg.RIU,
    source=source,
)


mesh = FibonacciMesh(
    sampling=1000,
    phi_offset=0 * ureg.degree,
    max_angle=10 * ureg.degree,
    gamma_offset=0 * ureg.degree,
)

# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)


res = experiment.binding._get_farfields(
    scatterer_set=experiment.scatterer.binding,
    source_set=experiment.source.binding,
    mesh=mesh,
)


print(res.shape)
