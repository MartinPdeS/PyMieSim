from PyMieSim.Source import GaussianBeam, PlaneWave
from PyMieSim.GLMT.Sphere import S1
from PyMieSim.GLMT.GaussianBeam import Anm, Anm_integrand
from PyMieSim.LMT.Sphere import Fields
import matplotlib.pyplot as plt
from mayavi import mlab
import matplotlib.pyplot as plt
from PyMieSim.Special import Pinm
import numpy as np
pi = np.pi
from ai import cs
from PyMieSim.Scatterer import Sphere

w = 0.6328e-6#1e-6;
k= 2*pi/w;
w0 = 2.0142649597710274e-06
s = 1/(k*w0)

print(s)
offset= [5e-6,3e-6,1e-16]
Offset = [5e-6*k,3e-6*k,1e-16*k]


Anm2 = Anm_integrand(n        = 20,
            m        = 2,
            sampling = 1000,
            k        = k,
            w0       = w0,
            s        = s,
            Offset   = Offset,
            offset   = offset,
            R0       = 57.89657280746359,
            xi       = 0.5404195002705842)


Anm15= Anm_integrand(n        = 20,
            m        = 15,
            sampling = 500,
            k        = k,
            w0       = w0,
            s        = s,
            Offset   = Offset,
            offset   = offset,
            R0       = 57.89657280746359,
            xi       = 0.5404195002705842)


fig = plt.figure(figsize=(8,5))
ax0 = fig.add_subplot(211)
ax1 = fig.add_subplot(212)
ax0.plot(Anm2.real, 'k', label=r"$\mathcal{Re}$ {$A_{nm}$}")
ax0.plot(Anm2.imag, 'r--',label=r"$\mathcal{Im}$ {$A_{nm}$}")

ax1.plot(Anm15.real, 'k', label=r"$\mathcal{Re}$ {$A_{nm}$}")
ax1.plot(Anm15.imag, 'r--',label=r"$\mathcal{Im}$ {$A_{nm}$}")
ax0.grid()
ax1.grid()
plt.show()







#-
