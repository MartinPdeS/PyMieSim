from PyMieSim.Source import GaussianBeam
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
pi = np.pi




Sampling=1000

beam = GaussianBeam(Wavelength = 0.6328e-6,
                    NA         = 0.142,
                    Offset     = [5e-6]*3)

mList = np.arange(0,50)
Anm20 = []; Anm180 = []
for m in mList:
    Anm20.append(np.abs( beam.Anm(20,m, Sampling=Sampling) ) )
    Anm180.append(np.abs( beam.Anm(180,m, Sampling=Sampling) ) )



dataStr0 = \
f"""\
n = 20
$x_0$ = {beam.offset[0]*1e6} $\mu$m
$y_0$ = {beam.offset[1]*1e6}  $\mu$m
$z_0$ = {beam.offset[2]*1e6} $\mu$m
$w_0$ = {beam.w0*1e6:.2f} $\mu$m\
"""

anchored_text0 = AnchoredText(dataStr0, loc=4)

dataStr1 = \
f"""\
n = 180
$x_0$ = {beam.offset[0]*1e6} $\mu$m
$y_0$ = {beam.offset[1]*1e6}  $\mu$m
$z_0$ = {beam.offset[2]*1e6} $\mu$m
$w_0$ = {beam.w0*1e6:.2f} $\mu$m\
"""

anchored_text1 = AnchoredText(dataStr1, loc=4)

fig =plt.figure(figsize=(8,5))
ax0 = fig.add_subplot(211)
ax1 = fig.add_subplot(212)

ax0.plot(mList, Anm20, 'C0s-')
ax1.plot(mList, Anm180, 'C1s-')

ax0.grid(); ax1.grid()
ax0.add_artist(anchored_text0)
ax1.add_artist(anchored_text1)

ax1.set_xlabel(r'$m$')
ax0.set_ylabel(r'$|A_{nm}|$')
ax1.set_ylabel(r'$|A_{nm}|$')
plt.show()



















#-
