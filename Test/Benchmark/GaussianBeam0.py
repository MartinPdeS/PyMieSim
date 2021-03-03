from PyMieSim.Source import GaussianBeam
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np


""" Results shoudl be the same as ref[2] figure 2. """

Theta = np.linspace(0,np.pi,20)

beam = GaussianBeam(Wavelength   = 0.6328e-6,
                    NA           = 0.2,
                    Polarization = 0,
                    Offset       = [5e-6,3e-6,0])


angle0, val0 = beam.Anm_integrand(20,2, Sampling=500)
angle1, val1 = beam.Anm_integrand(20,15, Sampling=500)

fig = plt.figure(figsize=(8,5))
ax0 = fig.add_subplot(211)
ax1 = fig.add_subplot(212)
ax0.plot(np.rad2deg(angle0), val0.real,'k', label=r"$\mathcal{Re}$ {$A_{nm}$}")
ax0.plot(np.rad2deg(angle0), val0.imag,'r--', label=r"$\mathcal{Im}$ {$A_{nm}$}")


ax1.plot(np.rad2deg(angle1), val1.real,'k', label=r"$\mathcal{Re}$ {$A_{nm}$}")
ax1.plot(np.rad2deg(angle1), val1.imag,'r--', label=r"$\mathcal{Im}$ {$A_{nm}$}")



dataStr0 = \
f"""\
n = 20; m = 2
$x_0$ = {beam.offset[0]*1e6} $\mu$m
$y_0$ = {beam.offset[1]*1e6}  $\mu$m
$z_0$ = {beam.offset[2]*1e6} $\mu$m
$w_0$ = {beam.w0*1e6:.2f} $\mu$m\
"""

anchored_text0 = AnchoredText(dataStr0, loc=4)


dataStr1 = \
f"""\
n = 20; m = 15
$x_0$ = {beam.offset[0]*1e6} $\mu$m
$y_0$ = {beam.offset[1]*1e6}  $\mu$m
$z_0$ = {beam.offset[2]*1e6} $\mu$m
$w_0$ = {beam.w0*1e6:.2f} $\mu$m\
"""

anchored_text1 = AnchoredText(dataStr1, loc=4)

ax0.grid()
ax1.grid()
ax0.legend(loc=3)
ax1.legend(loc=3)
ax0.add_artist(anchored_text0)
ax1.add_artist(anchored_text1)
ax0.set_ylabel(r'Integrand of $A_{nm}$')
ax1.set_ylabel(r'Integrand of $A_{nm}$')
ax1.set_xlabel(r'Angle $\theta$ [degree]')
plt.show()
