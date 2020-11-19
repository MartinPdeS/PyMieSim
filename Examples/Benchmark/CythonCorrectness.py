
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit


import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.cython.S1S2 import MieS1S2 as S1S2_CYTHON
from PyMieCoupling.functions.MieComputing import GetS1S2 as S1S2_PYTHON
AngleList = np.linspace(0,np.pi/2,200).tolist()


resPython =S1S2_PYTHON(1.4, 0.3, AngleList);

resCython = S1S2_CYTHON(1.4, 0.3, AngleList);

print(np.shape(resCython))
print(np.shape(resPython))
fig = plt.figure()
ax0 = fig.add_subplot(1,2,1)
ax1 = fig.add_subplot(1,2,2)
ax0.plot(resPython[1])
ax1.plot(resCython[1])
plt.show()
