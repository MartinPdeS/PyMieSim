
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
For different scatterer diameters.
_________________________________________________________
"""

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from src.classes.Detector import Detector
from src.classes.Scattering import Scatterer
from src.functions.couplings import PointFieldCoupling

npts=101

Detector = Detector(size       = 50e-6,
                    wavelength = 400e-9,
                    npts       = npts)

Detector.magnificate(magnification=1.5)

Detector.PlotFields()

DiameterList = np.linspace(100,9000,50) * 1e-9

Coupling = []

for Diameter in tqdm(DiameterList, total = len(DiameterList), desc ="Progress:"):

    Scat = Scatterer(diameter    = Diameter,
                     wavelength  = 400e-9,
                     index       = 1.4,
                     npts        = 101,
                     ThetaBound  = [-20,20],
                     ThetaOffset = 0,
                     PhiBound    = [-20,20],
                     PhiOffset   = 10)

    Coupling.append( PointFieldCoupling(Detector     = Detector,
                                            Source   = Scat.Field.Parallel,
                                            Mesh     = Scat.Meshes) )





fig = plt.figure(figsize=(15,5))
ax0 = fig.add_subplot(111)
ax0.plot(DiameterList*1e6, Coupling, 'C0', label=r'LP$_{01}$')
ax0.set_title('Mode coupling vs. Scatterer diameter')
ax0.legend()
ax0.set_xlabel(r'Scatterer diameter [$\mu m$]')
ax0.set_ylabel('Modal Coupling')
ax0.set_yscale('log')
ax0.grid()
plt.show()










# -
