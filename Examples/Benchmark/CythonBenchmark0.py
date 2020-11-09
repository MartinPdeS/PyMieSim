
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit

setup = """
from PyMieCoupling.functions.MieComputing import MieS1S2 as S1S2_PYTHON
from PyMieCoupling.S1S2 import MieS1S2 as S1S2_CYTHON"""


BenchPython = """
for i in range(10000):
    S1S2_PYTHON(1.4, 10, 0.5)"""

BenchCython = """
for i in range(10000):
    S1S2_CYTHON(1.4, 10, 0.5)"""


print('\nCYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchCython,
                    number = 1) )


print('='*50)

print('\nPYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchPython,
                    number = 1) )
