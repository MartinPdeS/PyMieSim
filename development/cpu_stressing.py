
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
import numpy as np
from PyMieSim.units import nanometer, degree, watt, AU, RIU

source = Gaussian(
    wavelength=400 * nanometer,
    polarization=[0] * degree,
    optical_power=1e-6 * watt,
    NA=0.2 * AU
)

scatterer = Sphere(
    diameter=np.linspace(300, 1000, 100) * nanometer,
    property=[1.2, 1.25] * RIU,
    medium_property=[1.0] * RIU,
    source=source
)

experiment = Setup(scatterer=scatterer, source=source)

dataframe = experiment.get('a1')

dataframe.plot_data(x='scatterer:diameter')
