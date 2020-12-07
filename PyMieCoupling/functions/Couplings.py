import numpy as np


""" Coupling Reference: Estimation of Coupling Efficiency of Optical Fiber by Far-Field Method """






def CenteredCoupling_Para(Detector, Scatterer):
    if Detector._coupling == "Intensity":
        Para = (Detector.Fourier * Scatterer.Field.Parallel).__abs__()**2
        Para = Para * np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian) )
        Para = Para.sum() * Detector.Meshes.dOmega.Radian


    if Detector._coupling == "Amplitude":
        Para = (Detector.Fourier * Scatterer.Field.Parallel)
        Para = Para * np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian) )
        Para = Para.sum() * Detector.Meshes.dOmega.Radian
        Para = Para.__abs__()**2

    return np.asscalar( Para )


def CenteredCoupling_Perp(Detector, Scatterer):
    if Detector._coupling == "Intensity":
        Perp = (Detector.Fourier * Scatterer.Field.Perpendicular ).__abs__()**2
        Perp = Perp * np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian) )
        Perp = Perp.sum() * Detector.Meshes.dOmega.Radian

    if Detector._coupling == "Amplitude":
        Perp = (Detector.Fourier * Scatterer.Field.Perpendicular)
        Perp = Perp * np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian) )
        Perp = Perp.sum()
        Perp = Perp * Detector.Meshes.dOmega.Radian
        Perp = Perp.__abs__()**2

    return np.asscalar( Perp )


def CenteredCoupling_Filtered(Detector, Scatterer):
    Perp = CenteredCoupling_Perp(Detector, Scatterer)
    Para = CenteredCoupling_Para(Detector, Scatterer)

    return  Perp * np.cos(Detector.FilterRad) + Para * np.sin(Detector.FilterRad)


def CenteredCoupling_NoFilter(Detector, Scatterer):
    Perp = CenteredCoupling_Perp(Detector, Scatterer)
    Para = CenteredCoupling_Para(Detector, Scatterer)

    return Perp + Para



def MeanCoupling_Para(Detector, Scatterer):
    Para = (Detector.Fourier * Scatterer.Field.Parallel).__abs__()**2
    Para = Para.sum()

    return np.asscalar( Para )


def MeanCoupling_Perp(Detector, Scatterer):
    Perp = (Detector.Fourier * Scatterer.Field.Perpendicular).__abs__()**2
    Perp = Perp.sum()

    return np.asscalar( Perp )


def GetFootprint(Detector, Scatterer):
    Perp = (Detector.Fourier * Scatterer.Field.Perpendicular).__abs__()**2
    Para = (Detector.Fourier * Scatterer.Field.Parallel).__abs__()**2

    return Perp, Para


def MeanCoupling_Filtered(Detector, Scatterer):
    Perp = MeanCoupling_Perp(Detector, Scatterer)
    Para = MeanCoupling_Para(Detector, Scatterer)

    return  Perp * np.cos(Detector.FilterRad)  + Para * np.sin(Detector.FilterRad)


def MeanCoupling_NoFilter(Detector, Scatterer):
    Perp = MeanCoupling_Perp(Detector, Scatterer)
    Para = MeanCoupling_Para(Detector, Scatterer)

    return Perp + Para



def Coupling(Scatterer, Detector, Polarization, Mode='Centered'):

    if Mode == 'Centered':
        return CenteredCoupling(Scatterer, Detector, Polarization)

    if Mode == 'Mean':
        return MeanCoupling(Scatterer, Detector, Polarization)



def CenteredCoupling(Scatterer, Detector, Polarization):

    if Polarization in ['Parallel']:
        return CenteredCoupling_Para(Detector, Scatterer)


    if Polarization in ['Perpendicular']:
        return CenteredCoupling_Perp(Detector, Scatterer)


    if Polarization in ['Filtered']:
        return CenteredCoupling_Filtered(Detector, Scatterer)


    if Polarization in ['NoFiltered']:
        return CenteredCoupling_NoFilter(Detector, Scatterer)


    if Polarization in ['all']:
        CouplingDict = {}

        coupling = CenteredCoupling_Para(Detector, Scatterer)
        CouplingDict['Parallel'] = coupling

        coupling = CenteredCoupling_Perp(Detector, Scatterer)
        CouplingDict['Perpendicular'] = coupling

        coupling = CenteredCoupling_Filtered(Detector, Scatterer)
        CouplingDict['Filtered'] = coupling

        coupling = CenteredCoupling_NoFilter(Detector, Scatterer)
        CouplingDict['NoFiltered'] = coupling

        return CouplingDict






def MeanCoupling(Scatterer, Detector, Polarization):

    if Polarization in ['Parallel']:
        return MeanCoupling_Para(Detector, Scatterer)

    if Polarization in ['Perpendicular']:
        return MeanCoupling_Perp(Detector, Scatterer)

    if Polarization in ['Filtered']:
        return MeanCoupling_Filtered(Detector, Scatterer)

    if Polarization in ['NoFiltered']:
        return MeanCoupling_NoFilter(Detector, Scatterer)


    if Polarization in ['all']:
        CouplingDict = {}

        coupling = MeanCoupling_Para(Detector, Scatterer)
        CouplingDict['Parallel'] = coupling

        coupling = MeanCoupling_Perp(Detector, Scatterer)
        CouplingDict['Perpendicular'] = coupling

        coupling = MeanCoupling_Filtered(Detector, Scatterer)
        CouplingDict['Filtered'] = coupling

        coupling = MeanCoupling_NoFilter(Detector, Scatterer)
        CouplingDict['NoFiltered'] = coupling

        return CouplingDict
