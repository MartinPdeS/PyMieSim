import numpy                    as np
from PyMieSim.Tools._Coupling   import *

def Coupling(Scatterer, Detector):
    EPhi, ETheta = Scatterer.uFarField(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian, 1)
    dOmega     = float(Detector.Mesh.dOmega.Radian)
    Omega      = float(Detector.Mesh.Omega.Radian)
    Filter     = Detector.Filter.Radian
    Scalar     = Detector.Scalar

    kwarg = {'ScalarField': Scalar.astype(np.complex),
             'EPhi'   :     EPhi.astype(np.complex),
             'ETheta' :     ETheta.astype(np.complex),
             'dOmega' :     dOmega}

    if Detector.CouplingMode[1] == 'Centered':
        if Detector.CouplingMode[0] == "Intensity":
            if Filter is None:
                return NoCoherentPointCoupling(**kwarg)
            else:
                return NoCoherentPointCouplingFilter(**kwarg, Filter = Filter)

        elif Detector.CouplingMode[0] == "Amplitude":
            if Filter is None:
                return CoherentPointCoupling(**kwarg)
            else:
                return CoherentPointCouplingFilter(**kwarg, Filter = Filter)


    elif Detector.CouplingMode[1] == 'Mean':
        if Detector.CouplingMode[0] == "Intensity":                                    # same thing as intensity point coupling
            if Filter is None:
                return NoCoherentMeanCoupling(**kwarg, Omega = Omega)
            else:
                return NoCoherentMeanCouplingFilter(**kwarg, Filter = Filter, Omega = Omega)



        elif Detector.CouplingMode[0] == "Amplitude":
            if Filter is None:
                return CoherentMeanCoupling(**kwarg, Omega = Omega)
            else:
                return CoherentMeanCouplingFilter(**kwarg, Filter = Filter, Omega = Omega)
