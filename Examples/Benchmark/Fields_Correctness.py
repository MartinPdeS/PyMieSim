
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""
import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.cython.S1S2 import GetS1S2 as S1S2_CYTHON
from PyMieCoupling.python.S1S2 import GetS1S2 as S1S2_PYTHON
from PyMieCoupling.cpp.S1S2 import GetFields as Fields_CPP
PhiList = np.linspace(0,np.pi/2,50)
ThetaList = np.linspace(0,np.pi/2,50)
PHI, THETA = np.meshgrid(ThetaList, PhiList)


S1S2_CYTHON = S1S2_CYTHON(1.4, 0.3, PhiList);
Parallel_CYTHON = np.outer(S1S2_CYTHON[0], np.sin(ThetaList))
Perpendicular_CYTHON = np.outer(S1S2_CYTHON[1], np.cos(ThetaList))



Parallel_CPP, Perpendicular_CPP = Fields_CPP(1.4, 0.3, ThetaList, PhiList, Polarization='None');
Parallel_CPP = Parallel_CPP.reshape([50,50])
Perpendicular_CPP = Perpendicular_CPP.reshape([50,50])

fig = plt.figure()
ax0 = fig.add_subplot(141)
ax1 = fig.add_subplot(142)
ax2 = fig.add_subplot(143)
ax3 = fig.add_subplot(144)

ax0.imshow(np.abs(Perpendicular_CYTHON))
ax1.imshow(np.abs(Perpendicular_CPP))
ax2.imshow(np.abs(Parallel_CYTHON))
ax3.imshow(np.abs(Parallel_CPP))

plt.show()


# 1
