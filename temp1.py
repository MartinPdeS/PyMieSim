from PyMieSim.Source import GaussianBeam, PlaneWave
from PyMieSim.GLMT.Sphere import S1
import matplotlib.pyplot as plt
from mayavi import mlab
import numpy as np
pi = np.pi
from ai import cs


Phi          = np.linspace(-pi/2,pi/2,100),
Theta        = np.linspace(-pi,pi,100),
PHI, THETA   = np.meshgrid(Phi, Theta)

s1 = S1(Index        = 1.4,
        Diameter     = 1e-6,
        Wavelength   = 1e-6,
        nMedium      = 1.0,
        Phi          = Phi,
        Theta        = Theta,
        Polarization = 0,
        E0           = 1,
        R            = 1,
        Lenght       = 100,
        BSC_TM       = {1:1, -1:-1},
        BSC_TE       = {1:1, -1:1})



x,y,z = cs.sp2cart(np.abs(s1).reshape([100,100]), PHI, THETA)
print(s1.shape)
mlab.mesh(x,y, z)
mlab.show()
