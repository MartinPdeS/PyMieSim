import numpy as np

from PyMieCoupling.utils import PlotUnstructureData

""" Coupling Reference: Estimation of Coupling Efficiency of Optical Fiber by Far-Field Method """






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


def GetFootprint(Detector, Scatterer, Num):

    MaxPhi = Detector.Mesh.Phi.Radian.max()

    Phi, Theta = Detector.StructuredSphericalMesh(Num = Num, MaxAngle=MaxPhi)

    Perp =  (Detector.StructuredFarField(Num=Num, SFactor=16) *\
    Scatterer.Perpendicular( Phi.flatten(), Theta.flatten() ).reshape(Theta.shape) )

    Para = (Detector.StructuredFarField(Num=Num, SFactor=16) *\
    Scatterer.Parallel( Phi.flatten(), Theta.flatten() ).reshape(Theta.shape) )#

    FourierPara = np.fft.ifft2(Para, s=[512*2, 512*2])

    FourierPara = np.fft.fftshift(FourierPara).__abs__()**2

    FourierPerp = np.fft.ifft2(Perp,  s=[512*2, 512*2])

    FourierPerp = np.fft.fftshift(FourierPerp).__abs__()**2

    return FourierPara + FourierPerp



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
