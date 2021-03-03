from PyMieSim.Source import GaussianBeam, PlaneWave
from PyMieSim.GLMT.Sphere import S1
from PyMieSim.LMT.Sphere import Fields
import matplotlib.pyplot as plt
from mayavi import mlab
import matplotlib.pyplot as plt
from PyMieSim.Special import Pinm
import numpy as np
pi = np.pi
from ai import cs



Num = 150
Phi          = np.linspace(0, pi, Num)
Theta        = np.linspace(0, 2*pi, Num)
PHI, THETA   = np.meshgrid(Phi, Theta)


beam = GaussianBeam(Wavelength=1e-6, NA = 0.6, Offset=[0e-6,0,0])
BSC = beam.GetBSC(Precision=1)
print(BSC)


s1, s2 = S1(Index    = 1.4,
        Diameter     = 2e-6,
        Wavelength   = 2e-6,
        nMedium      = 1.0,
        Phi          = Phi,
        Theta        = Theta,
        Polarization = 0,
        E0           = 1,
        R            = 1,
        BSC          = np.array(BSC.reset_index(level=[0,1]).values.tolist()),
        MaxOrder     = beam.NOrders
        )


temp = (np.abs(s1)**2 + np.abs(s2)**2)*1e20
print(np.abs(s1))

x,y,z = cs.sp2cart(temp.reshape(THETA.shape), PHI-pi/2, THETA)

mlab.mesh(x,y, z)
mlab.show()



















#-
