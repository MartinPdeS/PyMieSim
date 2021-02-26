from PyMieSim.Source import GaussianBeam
import matplotlib.pyplot as plt
import numpy as np

beam = GaussianBeam(Wavelength=0.632e-6,
                    NA           = 0.15,
                    Polarization = 0,
                    offset = [5e-6,3e-6,0e-6])

print(beam.w0)
val=[]
nList = np.arange(1,23)*10
for n in nList:
    val.append( beam.Anm(n,0) )


plt.semilogy(nList, np.abs(val),'k')

plt.grid()
plt.show()
