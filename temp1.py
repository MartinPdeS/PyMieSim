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



Num = 250
Phi          = np.linspace(0, pi, Num)
Theta        = np.linspace(0, 2*pi, Num)
PHI, THETA   = np.meshgrid(Phi, Theta)

MaxOrder     = 20

beam = PlaneWave(Wavelength=1e-6)
BSC = beam.GetBSC(MaxOrder=MaxOrder)



print(np.array(BSC.reset_index(level=[0,1]).values.tolist()))


s1,s2 = Fields(Index     = 1.4,
              Diameter     = 2e-6,
              Wavelength   = 1e-6,
              nMedium      = 1.0,
              Phi          = PHI.flatten(),
              Theta        = THETA.flatten(),
              Polarization = 0,
              E0           = 1,
              R            = 1,
              Lenght       = PHI.size
            )


temp = np.abs(s1)**2 + np.abs(s2)**2
x,y,z = cs.sp2cart(temp.reshape(THETA.shape), PHI-pi/2, THETA)

mlab.mesh(x,y, z)
mlab.show()



















#-
