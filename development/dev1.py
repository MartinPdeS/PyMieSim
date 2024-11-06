import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU

source = Gaussian(
    wavelength=np.linspace(400, 1000, 500) * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

scatterer = Sphere(
    diameter=[200] * nanometer,
    property=[4] * RIU,
    medium_property=1 * RIU,
    source=source
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get('a1', 'a2', 'a3', 'b1', 'b2', 'b3')

dataframe.plot_data(x="source:wavelength")
