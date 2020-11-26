
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.Optimization import LoopRIDiameter
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.DataFrame import Frame
from PyMieCoupling.classes.Optimizer import Simulator1 as Simulator

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)



Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


Detector0 = LPmode(Fiber         = Fiber,
                   Mode          = (0, 1),
                   Source        = LightSource,
                   Npts          = 51,
                   ThetaOffset   = 0,
                   PhiOffset     = 0,
                   NA            = 0.2)




Detector1 = Photodiode(NA                = 0.1,
                       Source            = LightSource,
                       Npts              = 101,
                       ThetaOffset       = 0,
                       PhiOffset         = 0)


Detector = Detector1

RIList = np.linspace(1.3, 1.5, 4).round(4)

DiameterList = np.linspace(500,1000,50).round(4) * 1e-9

def EvalFunc(x):

    Detector.PhiOffset = x

    Array = LoopRIDiameter(RIList       = RIList,
                           DiameterList = DiameterList,
                           Detector     = Detector,
                           Source       = LightSource,
                           QuietMode    = True,
                           Polarization = 'Filtered')    # can Parallel  -  Perpendicular  -  Filetered


    return Array.Cost('RI_RSD') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc)

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [10],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-6,
                  options  = {'maxiter': 50, 'rhobeg':20})
print(Result)

Detector.PhiOffset = Result.x


DF = Frame(RIList        = RIList,
           DiameterList  = DiameterList,
           Detector      = Detector,
           Source        = LightSource)

DF.Plot(y='Coupling')

#DF.Plot(y='STD')

plt.show()


print(Detector.Fourier)








# -
