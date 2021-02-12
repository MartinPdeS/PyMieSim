import matplotlib.pyplot as plt
import numpy as np

from PyMieSim.utils import PlotField
from PyMieSim.python.S1S2 import Field as PyField
from PyMieSim.cpp.Interface import GetFields as CppField

Phi = np.linspace(-np.pi/2,np.pi/2,50); Theta = np.linspace(-np.pi,np.pi,150)
PHI, THETA = np.meshgrid(Theta, Phi)


PyParallel, PyPerpendicular = PyField(1.4, 30, Phi-np.pi/2, Theta);

CppParallel, CppPerpendicular = CppField(1.4, 30, THETA.flatten(), PHI.flatten(), Polarization='None');

PlotField(Theta, Phi, CppParallel, CppPerpendicular)

PlotField(Theta, Phi, PyParallel, PyPerpendicular)


plt.show()



















    # -
