import GLMT
import numpy as np
from mayavi import mlab
from ai import cs


E0, E1, E2, R, Phi, Theta = GLMT.Sphere.ScatteredField(1e-6, 3, 1.4, 1, 1e-6, [1,2,3])

X, Y, Z = cs.sp2cart(E0, Phi, Theta)

mlab.mesh(X.real, Y.real, Z.real)
mlab.show()
