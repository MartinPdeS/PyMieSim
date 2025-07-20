
import numpy as np
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import nanometer, degree, RIU, volt, meter, radian

source = PlaneWave(
    wavelength=632.8 * nanometer,
    polarization=0 * degree,
    amplitude=1 * volt / meter,
)


scatterer = Sphere(
    diameter=200 * nanometer,
    property=1.5 * RIU,
    medium_property=1.0 * RIU,
    source=source,
)


phi = np.linspace(-180, 180) * degree


S1, S2 = scatterer.get_s1s2_array(phi)
import matplotlib.pyplot as plt

plt.figure()

plt.plot(phi, S1, label='S1')
plt.show()