"""
Goniometric Coupling vs S1 S2 Comparison
========================================

"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# PyMieSim imports
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment.scatterer import Sphere as ExperimentSphere
from PyMieSim.experiment.source import Gaussian as ExperimentGaussian
from PyMieSim.experiment import Setup
from PyMieSim.units import degree, watt, AU, RIU, micrometer

from PyMieSim.single.scatterer import Sphere as SingleSphere
from PyMieSim.single.source import Gaussian as SingleGaussian
from MPSPlots.styles import mps

# Setup parameters
scatterer_diameter = 0.3 * micrometer  # Diameter of the scatterer in meters
scatterer_index = 1.4 * RIU  # Refractive index of the scatterer
source_wavelength = 1.2 * micrometer  # Wavelength of the source in meters

# Experiment source and scatterer setup
source = ExperimentGaussian(
    wavelength=1.2 * micrometer,
    polarization=[0, 90] * degree,
    optical_power=1 * watt,
    NA=0.2 * AU
)

scatterer = ExperimentSphere(
    diameter=scatterer_diameter,
    index=scatterer_index,
    medium_index=1.0 * RIU,
    source=source
)

# Detector setup
detector = Photodiode(
    NA=[0.1] * AU,
    phi_offset=np.linspace(-180, 180, 100) * degree,
    gamma_offset=0.0 * degree,
    sampling=1000 * AU,
    polarization_filter=None
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
    polarization=90 * degree,
    optical_power=1 * watt,
    NA=0.2 * AU
)

single_scatterer = SingleSphere(
    diameter=scatterer_diameter,
    source=single_source,
    index=scatterer_index,
    medium_index=1.0 * RIU
)

s1s2 = single_scatterer.get_s1s2()
phi, s1, s2 = s1s2.phi, np.abs(s1s2.S1)**2, np.abs(s1s2.S2)**2
s1 /= s1.max()  # Normalize S1 data
s2 /= s2.max()  # Normalize S2 data


with plt.style.context(mps):
    figure, (ax0, ax1) = plt.subplots(1, 2, subplot_kw=dict(projection='polar'))



df = dataframe.unstack('polarization').pint.dequantify().reset_index().pint.quantify()
df['phi_offset'] /= (180 / np.pi)

df.plot(x='phi_offset', y=('coupling', 0 * degree), ax=ax0, linewidth=3, title='Polarization 90 degree')
df.plot(x='phi_offset', y=('coupling', 90 * degree), ax=ax1, linewidth=3, title='Polarization 90 degree')


ax0.plot(np.deg2rad(phi), s1, color='C1', linestyle='--', linewidth=1, label='Computed S1')

ax1.plot(np.deg2rad(phi), s2, color='C1', linestyle='--', linewidth=1, label='Computed S2')

plt.show()
