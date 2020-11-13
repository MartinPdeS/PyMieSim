
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import timeit

setup = """
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)"""


BenchPython = """
Scat = Scatterer(Diameter    = 10e-9,
                 Source      = LightSource,
                 Index       = 1.5,
                 Npts        = 201,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180],
                 CacheTrunk  = None,
                 cuda        = False)

Scat.S1S2                 """



print('\nCYTHON BENCHMARK')
print( timeit.timeit(setup = setup,
                    stmt = BenchPython,
                    number = 1000) )


print('='*50)
