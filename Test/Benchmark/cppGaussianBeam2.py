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

mList0 = np.arange(0,20)
mList1 = np.arange(0,50)
Anm20 = []; Anm180 = []
for m in mList0:
    print(m)
    Anm20.append( Anm(n      = 20,
                    m        = m,
                    k        = k,
                    w0       = w0,
                    s        = s,
                    Offset   = Offset,
                    offset   = offset,
                    R0       = 70.20990736659871,
                    xi       = 0.7853981633974484
                    ) )


for m in mList1:
    print(m)
    Anm180.append( Anm(n     = 180,
                    m        = m,
                    k        = k,
                    w0       = w0,
                    s        = s,
                    Offset   = Offset,
                    offset   = offset,
                    R0       = 70.20990736659871,
                    xi       = 0.7853981633974484
                    ) )

fig = plt.figure()
ax0 = fig.add_subplot(211)
ax1 = fig.add_subplot(212)
ax0.plot(mList0, np.abs(Anm20),'C0s-')
ax1.plot(mList1, np.abs(Anm180),'C0s-')
ax0.grid(); ax1.grid()
plt.show()





#-
