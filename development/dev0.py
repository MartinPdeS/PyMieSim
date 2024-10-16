from PyMieSim.experiment.source import PlaneWave
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, RIU
from PyOptik import Material
import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt


wavelength = np.linspace(200, 1200, 100)
diameter = 50
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

pymiesim_dataframe = experiment.get('Qsca', 'Qabs', 'Qext', add_units=False).reset_index()
pymiesim_dataframe.columns = pymiesim_dataframe.columns.get_level_values(0)
ax = pymiesim_dataframe.plot(y=['Qsca', 'Qabs', 'Qext'], x='source:wavelength', linewidth=4, figsize=(10, 6))

m = Material.gold.compute_refractive_index(wavelength * 1e-9)
_, qext_pymiescatt, qsca_pymiescatt, qabs_pymiescatt, *_  = PyMieScatt.MieQ_withWavelengthRange(m, diameter, n_medium, wavelength)


for name, value in zip(['Qsca', 'Qabs', 'Qext'], [qext_pymiescatt, qabs_pymiescatt, qsca_pymiescatt]):
    ax.plot(wavelength, value, ls='-', lw=0.5, marker='o', markersize=4, c='black', label=r"$Q_{name}^{PyMieScatt}$")

plt.show()






