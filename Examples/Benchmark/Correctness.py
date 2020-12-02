
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

AngleList = np.linspace(0,np.pi/2,101)#.tolist()

SizeParam, index = 0.2, 1.4

resPython = S1S2_PYTHON(index, SizeParam, AngleList);

resCython = S1S2_CYTHON(index, SizeParam, AngleList);

resCpp = S1S2_CPP(index, SizeParam, AngleList);
resCpp = S1S2_CPP(index, SizeParam, AngleList);
print(np.reshape( resCpp,[2,101] ))

fig = plt.figure(figsize=(10,5))
ax0 = fig.add_subplot(1,3,1); ax1 = fig.add_subplot(1,3,2); ax2 = fig.add_subplot(1,3,3);


ax0.set_title('PyMieScatt [PYTHON] result')
ax0.plot(np.real(resPython[0]),'C0', label='S1'); ax0.plot(np.real(resPython[1]), 'C1', label='S2')
ax0.grid()
ax0.legend()

ax1.set_title('CYTHON result')
ax1.plot(np.real(resCython[0]),'C0', label='S1'); ax1.plot(np.real(resCython[1]),'C1', label='S2');
ax1.grid()
ax1.legend()

ax2.set_title('CPP result')
ax2.plot(np.real(resCpp[0]),'C0', label='S1'); ax2.plot(np.real(resCpp[1]),'C1', label='S2');
ax2.grid()
ax2.legend()

plt.show()













#-
