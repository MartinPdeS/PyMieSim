import matplotlib.pyplot as plt
import numpy as np

from PyMieSim.utils import PlotField
from PyMieSim.LMT.python.Sphere import Fields as PyField
from PyMieSim.LMT.Sphere import GetFields as CppField


Phi = np.linspace(-np.pi/2,np.pi/2,50); Theta = np.linspace(-np.pi,np.pi,60)
PHI, THETA = np.meshgrid(Theta, Phi)


PyParallel, PyPerpendicular = PyField(1.4,
                                      1e-6,
                                      1e-6,
                                      1,
                                      Theta.flatten(),
                                      Phi.flatten(),
                                      0,
                                      1,
                                      1);

CppParallel, CppPerpendicular = CppField(1.4,
                                         1e-6,
                                         1e-6,
                                         1,
                                         THETA.flatten(),
                                         PHI.flatten(),
                                         0,
                                         1,
                                         1);

CppParallel = CppParallel.reshape([Phi.size, Theta.size])
CppPerpendicular = CppPerpendicular.reshape([Phi.size, Theta.size])


PlotField(Theta, Phi, CppParallel, CppPerpendicular)

PlotField(Theta, Phi, PyParallel, PyPerpendicular)


plt.show()



















    # -
