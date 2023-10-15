"""
CoreShell: Qsca, Qback vs core diameter
=======================================
"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.experiment import SourceSet, Setup, CoreShellSet
from PyMieSim.materials import BK7, Silver
from PyMieSim import measure

# %%
# Defining the ranging parameters for the scatterer distribution
# Here we look at core/shell scatterers and use constant shell diameter
# with variable core diameter
scatterer_set = CoreShellSet(
    core_diameter=numpy.geomspace(100e-09, 600e-9, 400),
    shell_width=800e-9,
    core_material=Silver,
    shell_material=BK7,
    n_medium=1
)

# %%
# Defining the source to be employed.
# The source is always a plane wave in the LMT framework.
# The amplitude is set to one per default.
source_set = SourceSet(
    wavelength=[800e-9, 900e-9, 1000e-9],
    polarization=0,
    amplitude=1
)

# %%
# Defining the experiment setup
# With the source and scatterers defined we set them together
# in an experiment.
experiment = Setup(
    scatterer_set=scatterer_set,
    source_set=source_set
)

# %%
# Measuring the properties
# We are interesting here in the back scattering efficientcy.
data = experiment.Get(measures=[measure.Qback])

# %%
# Plotting the results
figure = data.plot(
    y=measure.Qback,
    x=scatterer_set.core_diameter,
    y_scale='log'
)

_ = figure.show()
