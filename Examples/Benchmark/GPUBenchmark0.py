
"""
_________________________________________________________
Benchmark: comparison Cupy Vs. Numpy calling Stokes().
_________________________________________________________
"""

import timeit



SetupCPU = '''
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 200e-9,
                 Source      = LightSource,
                 Index       = 1.5,
                 Npts        = 401,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180],
                 GPU         = False,
                 CacheTrunk  = None) '''



SetupGPU = '''
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 200e-9,
                 Source      = LightSource,
                 Index       = 1.5,
                 Npts        = 401,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180],
                 GPU         = True,
                 CacheTrunk  = None) '''




Snippet = """
Scat.Stokes  """


print('\nCPU Benchmark:')
print(timeit.timeit(setup  = SetupCPU,
                    stmt   = Snippet,
                    number = 1000  ))
print('='*50)



print('\nGPU Benchmark:')
print(timeit.timeit(setup  = SetupGPU,
                    stmt   = Snippet,
                    number = 1000  ))
print('='*50)











# -
