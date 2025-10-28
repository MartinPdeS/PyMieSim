"""
Far-field scattering simulation for bulk scatterers
===================================================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
from TypedUnit import ureg
import numpy
import matplotlib.pyplot as plt
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.single.polarization import RightCircular
from PyOptik import Material
from PyMieSim.single.mesh import FibonacciMesh
from PyMieSim.utils import spherical_to_cartesian

# %%
# Defining the source to be employed.
source = Gaussian(
    wavelength=1200 * ureg.nanometer,
    polarization=RightCircular(),
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)
# %%
# Defining the ranging parameters for the scatterer distribution
scatterer = Sphere(
    diameter=[800, 1000, 1200] * ureg.nanometer,
    property=Material.BK7,
    medium_property=1 * ureg.RIU,
    source=source,
)

# %%
# Defining the detector to be employed.
mesh = FibonacciMesh(
    sampling=1400,
    max_angle=80 * ureg.degree,
    min_angle=0 * ureg.degree,
    phi_offset=0 * ureg.radian,
    gamma_offset=0 * ureg.radian,
)


# %%
# Defining the experiment setup
experiment = Setup(scatterer=scatterer, source=source)

# %%
# Measuring the properties
farfields = experiment.binding._get_farfields(
    scatterer_set=experiment.scatterer.set,
    source_set=experiment.source.set,
    mesh=mesh,
).squeeze()


spf = numpy.einsum("ijk->ik", abs(farfields) ** 2)


figure, axes = plt.subplots(
    ncols=3,
    nrows=1,
    figsize=(12, 4),
    subplot_kw={"projection": "3d"},
    sharex=True,
    sharey=True,
)


for idx, ax in enumerate(axes):

    x, y, z = spherical_to_cartesian(
        phi=mesh.spherical.phi,
        theta=mesh.spherical.theta,
        r=spf[idx],
    )

    ax.scatter(x, y, z)
    spf_min, spf_max = spf.min(), spf.max()
    ax.set(
        xlim=(-spf_max, spf_max),
        ylim=(-spf_max, spf_max),
        zlim=(-spf_max, spf_max),
    )


plt.show()
