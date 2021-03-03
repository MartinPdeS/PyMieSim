from PyMieSim.Source import GaussianBeam, PlaneWave
from PyMieSim.GLMT.Sphere import S1
from PyMieSim.LMT.Sphere import Fields
import matplotlib.pyplot as plt
from mayavi import mlab
import matplotlib.pyplot as plt
from PyMieSim.Special import Pinm
import numpy as np
pi = np.pi
from ai import cs
from PyMieSim.Scatterer import Sphere


Source = PlaneWave(Wavelength=1e-6)
scat = Sphere(Diameter=1e-6,Source = Source, Index=1.4)
print(scat.an(2))














#-
