
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit

setup = """
import numpy as np
from PyMieCoupling.cython.S1S2 import MieS1S2 as S1S2_CYTHON
from PyMieCoupling.python.S1S2 import MieS1S2 as S1S2_PYTHON
from PyMieCoupling.cpp.S1S2 import MieS1S2 as S1S2_CPP
AngleList = np.linspace(0,np.pi/2,200).tolist()"""


BenchPython = """S1S2_PYTHON(1.4, 0.3, AngleList);"""

BenchCython = """S1S2_CYTHON(1.4, 0.3, AngleList);"""

BenchCpp = """S1S2_CPP(1.4, 0.3, AngleList);"""


print('\nPYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchPython,
                    number = 100) )



print('='*50)


print('\nCYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchCython,
                    number = 100) )



print('='*50)


print('\nCPP BENCHMARK')
print( timeit.timeit(setup = setup,
                     stmt = BenchCpp,
                     number = 100) )









# 1
