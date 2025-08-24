"""
Cylinder Scatterer Bohren-Huffman figure 8.10
=============================================

"""

# %%
# Importing the dependencies: numpy, matplotlib, PyMieSim
import numpy
import matplotlib.pyplot as plt
from TypedUnit import ureg
from MPSPlots.styles import mps

from PyMieSim.directories import validation_data_path
from PyMieSim.single.source import Gaussian
from PyMieSim.single.scatterer import Cylinder

theoretical = numpy.genfromtxt(f"{validation_data_path}/bohren_huffman/figure_810.csv", delimiter=',')

x = theoretical[:, 0]
y = theoretical[:, 1]

source = Gaussian(
    wavelength=470 * ureg.nanometer,
    polarization=90 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.1 * ureg.AU,
)

scatterer = Cylinder(
    diameter=3000 * ureg.nanometer,
    source=source,
    property=(1.0 + 0.07j) * ureg.RIU,
    medium_property=1.0 * ureg.RIU
)

S1S2 = scatterer.get_s1s2(sampling=800)
data = (numpy.abs(S1S2.S1)**2 + numpy.abs(S1S2.S2)**2) * (0.5 / (numpy.pi * source.wavenumber.to_base_units()))**(1 / 4)

with plt.style.context(mps):
    figure, ax = plt.subplots(1, 1)

ax.plot(S1S2.phi, data, 'C1-', linewidth=3, label='PyMieSim')
ax.plot(x, y, 'k--', linewidth=1, label='B&H [8.10]')

ax.set(
    xlabel='scatterertering angle [degree]',
    ylabel='Phase function',
    yscale='log',
)

ax.legend()
plt.show()


# -
