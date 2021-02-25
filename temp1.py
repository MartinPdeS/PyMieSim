from PyMieSim.Source import GaussianBeam
import matplotlib.pyplot as plt
import numpy as np

beam = GaussianBeam(Wavelength=0.632e-6,
                    NA           = 0.8,
                    Polarization = 0)


x = 1e-6; y=1e-6; z=0
val= []
nList = np.arange(1,23) * 10
print(nList)
for n in nList:
    print(n)
    val.append(np.abs(beam.Anm(n,0)))

plt.semilogy(nList, val)
plt.grid()
plt.show()
