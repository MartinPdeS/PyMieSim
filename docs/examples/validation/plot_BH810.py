"""
Bohren-Huffman (figure~8.10)
============================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim
import numpy
import matplotlib.pyplot as plt

from PyMieSim.directories import validation_data_path
from PyMieSim.single.source import Gaussian
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.units import degree, watt, AU, RIU, nanometer

theoretical = numpy.genfromtxt(f"{validation_data_path}/Figure810BH.csv", delimiter=',')

x = theoretical[:, 0]
y = theoretical[:, 1]

source = Gaussian(
    wavelength=470 * nanometer,
    polarization=90 * degree,
    optical_power=1e-3 * watt,
    NA=0.1 * AU,
)

scatterer = Cylinder(
    diameter=3000 * nanometer,
    source=source,
    index=(1.0 + 0.07j) * RIU,
    medium_index=1.0 * RIU
)

S1S2 = scatterer.get_s1s2(sampling=800)
data = (numpy.abs(S1S2.S1)**2 + numpy.abs(S1S2.S2)**2) * (0.5 / (numpy.pi * source.wavenumber.to_base_units()))**(1 / 4)

plt.figure(figsize=(8, 4))
plt.plot(S1S2.phi, data, 'C1-', linewidth=3, label='PyMieSim')

plt.plot(x, y, 'k--', linewidth=1, label='B&H [8.10]')

plt.xlabel('scatterertering angle [degree]')
plt.ylabel('Phase function')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# -
