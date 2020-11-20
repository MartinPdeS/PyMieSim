
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit


import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.cython.S1S2 import MieS1S2 as S1S2_CYTHON
from PyMieCoupling.python.S1S2 import MieS1S2 as S1S2_PYTHON
from PyMieCoupling.cpp.S1S2 import MieS1S2 as S1S2_CPP

AngleList = np.linspace(0,np.pi/2,20).tolist()

resPython = S1S2_PYTHON(1.4, 0.3, AngleList);

resCython = S1S2_CYTHON(1.4, 0.3, AngleList);

resCpp = S1S2_CPP(1.4, 0.3, AngleList);


fig = plt.figure(figsize=(10,5))
ax0 = fig.add_subplot(1,3,1); ax1 = fig.add_subplot(1,3,2); ax2 = fig.add_subplot(1,3,3);

ax0.set_title('CYTHON result')
ax0.plot(np.abs(resPython[0]),'C0', label='S1'); ax0.plot(np.abs(resPython[1]), 'C1', label='S2')
ax0.grid()
ax0.legend()

ax1.set_title('PyMieScatt [PYTHON] result')
ax1.plot(np.abs(resCython[0]),'C0', label='S1'); ax1.plot(np.abs(resCython[1]),'C1', label='S2');
ax1.grid()
ax1.legend()

ax2.set_title('CPP result')
ax2.plot(np.abs(resCpp[0]),'C0', label='S1'); ax2.plot(np.abs(resCpp[1]),'C1', label='S2');
ax2.grid()
ax2.legend()

plt.show()













#-
