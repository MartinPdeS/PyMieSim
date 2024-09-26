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

from PyMieSim.experiment import measure
from PyMieSim.single.scatterer import Sphere as SingleSphere
from PyMieSim.single.source import Gaussian as SingleGaussian

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
experiment = Setup(
    scatterer=scatterer,
    source=source,
    detector=detector
)

# Gather data
coupling_data = experiment.get(measure.coupling, export_as='numpy').squeeze()
coupling_data /= coupling_data.max()  # Normalize data

# Single scatterer simulation for S1 and S2
single_source = SingleGaussian(
    wavelength=source_wavelength.to_base_units().magnitude,
    polarization=90,
    optical_power=1,
    NA=0.2
)

single_scatterer = SingleSphere(
    diameter=scatterer_diameter.to_base_units().magnitude,
    source=single_source,
    index=scatterer_index.to_base_units().magnitude,
    medium_index=1.0
)

s1s2 = single_scatterer.get_s1s2()
phi, s1, s2 = s1s2.phi, np.abs(s1s2.S1)**2, np.abs(s1s2.S2)**2
s1 /= s1.max()  # Normalize S1 data
s2 /= s2.max()  # Normalize S2 data


# Plotting the results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), subplot_kw={'projection': 'polar'})
ax1.plot(np.deg2rad(phi), s1, 'b-', linewidth=3, label='Computed S1')

ax1.plot(detector.phi_offset, coupling_data[0], 'k--', linewidth=1, label='Coupling for polarization: 0')

ax2.plot(np.deg2rad(phi), s2, 'r-', linewidth=3, label='Computed S2')
ax2.plot(detector.phi_offset, coupling_data[1], 'k--', linewidth=1, label='Coupling for polarization: 90')

ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
plt.tight_layout()
plt.show()
