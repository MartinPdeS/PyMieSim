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

theoretical = numpy.genfromtxt(f"{validation_data_path}/Figure810BH.csv", delimiter=',')

x = theoretical[:, 0]
y = theoretical[:, 1]

source = Gaussian(
    wavelength=470e-9,
    polarization=90,
    optical_power=1e-3,
    NA=0.1,
)

scatterer = Cylinder(
    diameter=3000e-9,
    source=source,
    index=1.0 + 0.07j,
    medium_index=1.0
)

S1S2 = scatterer.get_s1s2(sampling=800)
Data = (numpy.abs(S1S2.S1)**2 + numpy.abs(S1S2.S2)**2) * (0.5 / (numpy.pi * source.wavenumber))**(1 / 4)

plt.figure(figsize=(8, 4))
plt.plot(S1S2.phi, Data, 'C1-', linewidth=3, label='PyMieSim')

plt.plot(x, y, 'k--', linewidth=1, label='B&H [8.10]')

plt.xlabel('scatterertering angle [degree]')
plt.ylabel('Phase function')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# -
