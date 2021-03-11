from PyMieSim.Source import GaussianBeam
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
pi = np.pi




Sampling=1000

beam = GaussianBeam(Wavelength = 0.6328e-6,
                    NA         = 0.14,
                    Offset     = [5e-6]*3)


print(beam.w0, beam.Wavelength, beam.k, beam.R0, beam.xi)
nList = np.arange(1,150)
Anm = [];
for n in nList:
    print(n)
    Anm.append(np.abs( beam.Anm(n,0, Sampling=Sampling) ) )




dataStr = \
f"""\
n = 20
$x_0$ = {beam.offset[0]*1e6} $\mu$m
$y_0$ = {beam.offset[1]*1e6}  $\mu$m
$z_0$ = {beam.offset[2]*1e6} $\mu$m
$w_0$ = {beam.w0*1e6:.2f} $\mu$m\
"""

anchored_text = AnchoredText(dataStr, loc=4)

fig = plt.figure(figsize=(8,5))
ax0 = fig.add_subplot(111)

ax0.semilogy(nList, Anm, 'C0s-')

ax0.grid()
ax0.add_artist(anchored_text)


ax0.set_xlabel(r'$m$')
ax0.set_ylabel(r'$|A_{nm}|$')
plt.show()

















#-
