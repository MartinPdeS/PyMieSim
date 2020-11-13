
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit

setup = """
import numpy as np
from PyMieCoupling.functions.MieComputing import MieS1S2 as S1S2_PYTHON
from PyMieCoupling.S1S2 import MieS1S2 as S1S2_CYTHON
AngleList = np.linspace(0,np.pi/2,200)"""


BenchPython = """
for i in range(10):
    S1S2_PYTHON(1.4, 0.3, AngleList)"""

BenchCython = """
for i in range(10):
    a,b = S1S2_CYTHON(1.4, 0.3, AngleList)"""


print('\nCYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchCython,
                    number = 1) )


print('='*50)

print('\nPYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchPython,
                    number = 1) )
