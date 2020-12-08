
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.functions.converts import Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Detector import fiber, LPmode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Scattering import Scatterer


LightSource = Source(Wavelength   = 940e-9,
                     Polarization = 0)

npts=101

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = LPmode(Fiber         = Fiber,
              Mode          = (1, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Filter        = 0,
              NA            = 0.5)

LP01 = LPmode(Fiber         = Fiber,
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Filter        = 0,
              NA            = 0.5)


Scat = Scatterer(Diameter    = 100e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Meshes      = LP01.Meshes)

LP01.Fourier.Plot()
plt.show()
"""
Footprint = Scat.Footprint(LP01)
x = Angle2Direct(Scat.Meshes.Phi.Vector.Degree, LightSource.k)*1e6
y = Angle2Direct(Scat.Meshes.Theta.Vector.Degree, LightSource.k)*1e6

fig = plt.figure()
ax = fig.add_subplot(111)
ax.pcolormesh(x, y, Footprint, cmap='gray')
ax.set_aspect('equal')
ax.set_xlabel(r'X-Distance [$\mu m$]')
ax.set_ylabel(r'Y-Distance [$\mu m$]')

plt.show()

"""
# -
