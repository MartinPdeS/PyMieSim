"""
Goniometric Coupling vs S1 S2 Comparison
========================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
from TypedUnit import ureg

# PyMieSim imports
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere as ExperimentSphere
from PyMieSim.experiment.source import Gaussian as ExperimentGaussian
from PyMieSim.experiment import Setup

from PyMieSim.single.scatterer import Sphere as SingleSphere
from PyMieSim.single.source import Gaussian as SingleGaussian
from MPSPlots.styles import mps

# Setup parameters
scatterer_diameter = 0.3 * ureg.micrometer  # Diameter of the scatterer in meters
scatterer_index = 1.4 * ureg.RIU  # Refractive index of the scatterer
source_wavelength = 1.2 * ureg.micrometer  # Wavelength of the source in meters

# Experiment source and scatterer setup
source = ExperimentGaussian(
    wavelength=1.2 * ureg.micrometer,
    polarization=[0, 90] * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.2 * ureg.AU
)

scatterer = ExperimentSphere(
    diameter=scatterer_diameter,
    property=scatterer_index,
    medium_property=1.0 * ureg.RIU,
    source=source
)

# Detector setup
detector = Photodiode(
    NA=[0.1] * ureg.AU,
    phi_offset=np.linspace(-180, 180, 100) * ureg.degree,
    gamma_offset=0.0 * ureg.degree,
    sampling=1000 * ureg.AU,
)

# Configure experiment
experiment = Setup(scatterer=scatterer, source=source, detector=detector)

# Gather data
dataframe = experiment.get('coupling', drop_unique_level=True)
# dataframe.index = /= 180 / np.pi
dataframe['coupling'] /= dataframe['coupling'].max()  # Normalize data

# Single scatterer simulation for S1 and S2
single_source = SingleGaussian(
    wavelength=source_wavelength,
    polarization=90 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.2 * ureg.AU
)

single_scatterer = SingleSphere(
    diameter=scatterer_diameter,
    source=single_source,
    property=scatterer_index,
    medium_property=1.0 * ureg.RIU
)

s1s2 = single_scatterer.get_s1s2()
phi, s1, s2 = s1s2.phi, np.abs(s1s2.S1)**2, np.abs(s1s2.S2)**2
s1 /= s1.max()  # Normalize S1 data
s2 /= s2.max()  # Normalize S2 data


with plt.style.context(mps):
    figure, ax0 = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))


df = dataframe.unstack('source:polarization').pint.dequantify().reset_index().pint.quantify()
df['detector:phi_offset'] /= (180 / np.pi)

df.plot(x='detector:phi_offset', y='coupling', ax=ax0, linewidth=3, title='Polarization 90 degree')

ax0.plot(np.deg2rad(phi), s1, color='white', linestyle='--', linewidth=1, label='Computed S1')

ax0.plot(np.deg2rad(phi), s2, color='cyan', linestyle='--', linewidth=1, label='Computed S2')

ax0.grid()

plt.show()
