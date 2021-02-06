import numpy as np
from ai import cs


from PyMieCoupling.utils import PlotUnstructureData
from PyMieCoupling.functions.converts import Angle2Direct
from PyMieCoupling.classes.Representations import Footprint

""" Coupling Reference: Estimation of Coupling Efficiency of Optical Fiber by Far-Field Method """



def Direct2Angle(X, Y, MaxAngle):
    MaxZ = np.max(50) / np.tan(MaxAngle)

    Gamma = np.arctan(np.abs(X / MaxZ) )

    Phi = np.arctan(np.abs(X / MaxZ) )

    return Phi, Gamma


def AngleUnit2DirectUnit(Angle, k):
    FourierSpace = np.sin(Angle) * k / (2 * np.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = np.fft.fftshift( np.fft.fftfreq( Angle.shape[0], d = fourier_unit ) )

    return DirectSpace



def _CenteredCoupling(Detector, Scatterer):
    if Detector._CouplingMode[0] == "Intensity":
        Para = (Detector.Scalar * Scatterer.Parallel(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian)).__abs__()**2
        Para = Para * Detector.Mesh.SinMesh
        Para = Para.sum() * Detector.Mesh.dOmega.Radian

        Perp = (Detector.Scalar * Scatterer.Perpendicular(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian) ).__abs__()**2
        Perp = Perp * Detector.Mesh.SinMesh
        Perp = Perp.sum() * Detector.Mesh.dOmega.Radian

    if Detector._CouplingMode[0] == "Amplitude":
        Para = (Detector.Scalar * Scatterer.Parallel(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian))
        Para = Para * Detector.Mesh.SinMesh
        Para = Para.sum() * Detector.Mesh.dOmega.Radian
        Para = Para.__abs__()**2

        Perp = (Detector.Scalar * Scatterer.Perpendicular(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian))
        Perp = Perp * Detector.Mesh.SinMesh
        Perp = Perp.sum() * Detector.Mesh.dOmega.Radian
        Perp = Perp.__abs__()**2

    return np.asscalar( Para ), np.asscalar( Perp )



def MeanCoupling_Para(Detector, Scatterer):
    Para = (Detector.Scalar * Scatterer.Parallel(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian) * Detector.Mesh.SinMesh).__abs__()**2
    Para = Para.sum()

    return np.asscalar( Para )


def MeanCoupling_Perp(Detector, Scatterer):
    Perp = (Detector.Scalar * Scatterer.Perpendicular(Detector.Mesh.Phi.Radian, Detector.Mesh.Theta.Radian) * Detector.Mesh.SinMesh).__abs__()**2
    Perp = Perp.sum()

    return np.asscalar( Perp )

def Direct2spherical(X, Y, MaxAngle):
    Z = 50 / np.tan(MaxAngle)

    _, Phi, Theta = cs.cart2sp(X, Y, X*0+Z)

    return Phi, Theta

def Direct2Angle(X, Y, MaxAngle):
    MaxZ = np.max(X) / np.cos(MaxAngle)


def GetFootprint(Detector, Scatterer, Num):

    Footprin = Footprint(Scatterer = Scatterer, Detector = Detector, Num=200)

    return Footprin



def Coupling(Scatterer, Detector):
    if Detector._CouplingMode[1] == 'Centered':
        return CenteredCoupling(Scatterer, Detector)

    if Detector._CouplingMode[1] == 'Mean':
        return MeanCoupling(Scatterer, Detector)



def CenteredCoupling(Scatterer, Detector):

    if Detector._Filter.Radian != 'None':
        Para, Perp = _CenteredCoupling(Detector, Scatterer)
        Para *= np.sin(Detector._Filter.Radian)**2
        Perp *= np.cos(Detector._Filter.Radian)**2

    else:
        Para, Perp = _CenteredCoupling(Detector, Scatterer)

    return Perp + Para




def MeanCoupling(Scatterer, Detector):
    if Detector._Filter.Radian != 'None':
        Perp = MeanCoupling_Perp(Detector, Scatterer) * np.cos(Detector.Filter.Radian)**2
        Para = MeanCoupling_Para(Detector, Scatterer) * np.sin(Detector.Filter.Radian)**2

    else:
        Perp = MeanCoupling_Perp(Detector, Scatterer)
        Para = MeanCoupling_Para(Detector, Scatterer)

    return Perp + Para
