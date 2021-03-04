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

w = 6.328e-07
k= 9929180.321080256

w0 = 2.8775213711014675e-06
s = 1/(k*w0)

print(s)
offset= [5.e-06, 5.e-06, 5.e-06]
Offset = [49.64590161, 49.64590161, 49.64590161]

nList = np.arange(1,100)
Anm2 = []
for n in nList:
    print(n)
    Anm2.append( Anm(n       = n,
                    m        = 0,
                    k        = k,
                    w0       = w0,
                    s        = s,
                    Offset   = Offset,
                    offset   = offset,
                    R0       = 70.20990736659871,
                    xi       = 0.7853981633974484
                    ) )

print(np.abs(Anm2))
plt.figure()
plt.semilogy(nList, np.abs(Anm2))
plt.grid()
plt.show()





#-
