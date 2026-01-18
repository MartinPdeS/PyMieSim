import numpy as np
import matplotlib.pyplot as plt

from PyMieSim.units import ureg
from PyMieSim import single

ns = 1.60
NA = 0.95

source = single.source.Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=30 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

scatterer = single.scatterer.Sphere(
    diameter=500 * ureg.nanometer,
    source=source,
    property=(1.8 + 0.02j) * ureg.RIU,
    medium_property=ns * ureg.RIU,
)

scatterer.get_stokes(distance=2 * ureg.meter, sampling=80)

nd_values = np.linspace(1.0, 1.6, 13)
couplings = []

for nd in nd_values:
    detector = single.detector.Photodiode(
        NA=float(NA) * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=6000,
        medium_refractive_index=float(nd) * ureg.RIU,
    )
    couplings.append(float(detector.get_coupling(scatterer=scatterer).to("watt").magnitude))

plt.figure()
plt.plot(nd_values, couplings, marker="o")
plt.xlabel("Detector medium refractive index nd")
plt.ylabel("Coupling (W)")
plt.title(f"Edge test at fixed NA={NA}: coupling must vary with nd if interface mapping is implemented")
plt.grid(True)
plt.show()

print("nd:", nd_values.tolist())
print("coupling:", couplings)
