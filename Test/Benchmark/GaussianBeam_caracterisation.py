from PyMieSim.Source import GaussianBeam
import matplotlib.pyplot as plt
import numpy as np


""" Results shoudl be the same as ref[2] figure 2. """

beam = GaussianBeam(Wavelength=0.6328e-6,
                    NA           = 0.21,
                    Polarization = 0,
                    offset       = [0e-6,0e-6,0])


angle0, val0 = beam.Anm_integrand(20,2)
angle1, val1 = beam.Anm_integrand(20,15)

fig = plt.figure(figsize=(8,5))
ax0 = fig.add_subplot(211)
ax1 = fig.add_subplot(212)
ax0.plot(np.rad2deg(angle0), val0.real,'k', label=r"$\mathcal{Real}$ {$A_{nm}$}")
ax0.plot(np.rad2deg(angle0), val0.imag,'r--', label=r"$\mathcal{Imag}$ {$A_{nm}$}")


ax1.plot(np.rad2deg(angle1), val1.real,'k', label=r"$\mathcal{Real}$ {$A_{nm}$}")
ax1.plot(np.rad2deg(angle1), val1.imag,'r--', label=r"$\mathcal{Imag}$ {$A_{nm}$}")

ax0.grid()
ax1.grid()
ax0.legend()
ax1.legend()
ax0.set_ylabel(r'Integrand of $A_{nm}$')
ax1.set_ylabel(r'Integrand of $A_{nm}$')
ax1.set_xlabel(r'Angle $\theta$ [degree]')
plt.show()
