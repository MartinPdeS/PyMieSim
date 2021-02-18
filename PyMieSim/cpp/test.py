import GLMT
import numpy as np
from mayavi import mlab
from ai import cs
from numpy import array, arange, random
from tvtk.api import tvtk
from mayavi.scripts import mayavi2
num = 50

E0, E1, E2, R, Phi, Theta = GLMT.Sphere.ScatteredField(0.05e-6,
                                                       3,
                                                       1.4,
                                                       1.,
                                                       1e-6,
                                                       num)
print('PYTHON FINISH')


Theta = np.asarray(Theta).reshape([num,num])
Phi   = np.asarray(Phi).reshape([num,num])
R     = np.asarray(R).reshape([num,num])
E0    = np.asarray(E0).reshape([num,num])
E1    = np.asarray(E1).reshape([num,num])
E2    = np.asarray(E2).reshape([num,num])


E =  np.sqrt(np.abs(E0)**2 + np.abs(E1)**2 + np.abs(E2)**2)

print(E)

Theta = Theta

X, Y, Z = cs.sp2cart(E, Phi, Theta)

mlab.mesh(X, Y, Z)
#mlab.points3d(X, Y, Z, np.abs(X)**2 + np.abs(E1)**2 + np.abs(E2)**2)
mlab.show()
