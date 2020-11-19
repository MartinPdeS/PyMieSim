
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit

setup = """
import numpy as np
from PyMieCoupling.cython.S1S2 import MieS1S2 as S1S2_CYTHON
from PyMieCoupling.functions.MieComputing import GetS1S2 as S1S2_PYTHON
AngleList = np.linspace(0,np.pi/2,200).tolist()"""


BenchPython = """S1S2_PYTHON(1.4, 0.3, AngleList);"""

BenchCython = """S1S2_CYTHON(1.4, 0.3, AngleList);"""


print('\nCYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchCython,
                    number = 1) )


print('='*50)

print('\nPYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchPython,
                    number = 1) )
