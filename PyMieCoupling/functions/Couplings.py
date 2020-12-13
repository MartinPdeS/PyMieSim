import numpy as np
import matplotlib.pyplot as plt

""" Coupling Reference: Estimation of Coupling Efficiency of Optical Fiber by Far-Field Method """






def CenteredCoupling_Para(Detector, Scatterer):
    if Detector._coupling == "Intensity":
        Para = (Detector.FarField.Spherical * Scatterer.Field.Parallel).__abs__()**2
        Para = Para * np.sin(Detector.FarField.Meshes.Phi.Mesh.Radian.T).__abs__()
        Para = Para.sum() * Detector.FarField.Meshes.dOmega.Radian



    if Detector._coupling == "Amplitude":
        Para = (Detector.FarField.Spherical * Scatterer.Field.Parallel)
        Para = Para * np.sin(Detector.FarField.Meshes.Phi.Mesh.Radian.T).__abs__()
        Para = Para.sum() * Detector.FarField.Meshes.dOmega.Radian
        Para = Para.__abs__()**2

    return np.asscalar( Para )


def CenteredCoupling_Perp(Detector, Scatterer):
    if Detector._coupling == "Intensity":
        Perp = (Detector.FarField.Spherical * Scatterer.Field.Perpendicular ).__abs__()**2
        Perp = Perp * np.sin(Detector.FarField.Meshes.Phi.Mesh.Radian.T).__abs__()
        Perp = Perp.sum() * Detector.FarField.Meshes.dOmega.Radian

    if Detector._coupling == "Amplitude":

        Perp = (Detector.FarField.Spherical * Scatterer.Field.Perpendicular)
        Perp = Perp * np.sin(Detector.FarField.Meshes.Phi.Mesh.Radian.T).__abs__()
        Perp = Perp.sum()
        Perp = Perp * Detector.FarField.Meshes.dOmega.Radian
        Perp = Perp.__abs__()**2

    return np.asscalar( Perp )



def MeanCoupling_Para(Detector, Scatterer):
    Para = (Detector.FarField.Spherical * Scatterer.Field.Parallel * np.sin(Detector.FarField.Meshes.Phi.Mesh.Radian.T)).__abs__()**2
    Para = Para.sum()

    return np.asscalar( Para )


def MeanCoupling_Perp(Detector, Scatterer):
    Perp = (Detector.FarField.Spherical * Scatterer.Field.Perpendicular * np.sin(Detector.FarField.Meshes.Phi.Mesh.Radian.T)).__abs__()**2
    Perp = Perp.sum()

    return np.asscalar( Perp )


def GetFootprint(Detector, Scatterer):
    Perp = (Detector.FarField.Spherical * Scatterer.Field.Perpendicular).__abs__()**2
    Para = (Detector.FarField.Spherical * Scatterer.Field.Parallel).__abs__()**2

    return Perp + Para



def Coupling(Scatterer, Detector, Mode='Centered'):

    if Mode == 'Centered':
        return CenteredCoupling(Scatterer, Detector)

    if Mode == 'Mean':
        return MeanCoupling(Scatterer, Detector)



def CenteredCoupling(Scatterer, Detector):

    if Detector._Filter.Radian != 'None':
        Perp = CenteredCoupling_Perp(Detector, Scatterer) * (np.cos(Detector._Filter.Radian)**2)
        Para = CenteredCoupling_Para(Detector, Scatterer) * (np.sin(Detector._Filter.Radian)**2)

    else:
        Perp = CenteredCoupling_Perp(Detector, Scatterer)
        Para = CenteredCoupling_Para(Detector, Scatterer)


    return Perp + Para




def MeanCoupling(Scatterer, Detector):
    if Detector._Filter.Radian != 'None':
        Perp = MeanCoupling_Perp(Detector, Scatterer) * np.cos(Detector.Filter.Radian)**2
        Para = MeanCoupling_Para(Detector, Scatterer) * np.sin(Detector.Filter.Radian)**2
    else:
        Perp = MeanCoupling_Perp(Detector, Scatterer)
        Para = MeanCoupling_Para(Detector, Scatterer)

    return Perp + Para
