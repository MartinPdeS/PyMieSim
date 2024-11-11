from PyMieSim.experiment.source import PlaneWave
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, RIU
from PyOptik import Material
import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt


wavelength = np.linspace(200, 1200, 10)

diameter = [50, 100, 150]
n_medium = 1.5
material = Material.gold

#%% PyMieSim
source = PlaneWave(
    wavelength=wavelength * nanometer,
    polarization=0 * degree,
    amplitude=0 * watt,
)

scatterer = Sphere(
    diameter=diameter * nanometer,
    property=material,
    medium_property=n_medium * RIU,
    source=source
)

experiment = Setup(scatterer=scatterer, source=source)


dataframe = experiment.get('Qsca', add_units=False)

data = dataframe.values


data = data.reshape([10, 3])

print(data.shape)



# /
# print(res)


