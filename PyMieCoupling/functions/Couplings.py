import numpy as np



def Coupling_Para(Detector, Scatterer):

    Para = Detector.Fourier *\
           (Scatterer.Field.Parallel).__abs__() *\
           (np.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()

    return np.asscalar( ( Para.sum().__abs__() * Detector.Meshes.dOmega.Radian ) **2 )


def Coupling_Perp(Detector, Scatterer):

    Perp = Detector.Fourier *\
           (Scatterer.Field.Perpendicular).__abs__() *\
           (np.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()

    return np.asscalar( ( Perp.sum().__abs__() * Detector.Meshes.dOmega.Radian )**2 )


def Coupling_Filtered(Detector, Scatterer):

    Perp = Detector.Fourier *\
           (Scatterer.Field.Perpendicular).__abs__() *\
           (np.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()


    Para = Detector.Fourier *\
           (Scatterer.Field.Parallel).__abs__() *\
           (np.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()

    PerpFiltre = ( ( np.cos(Detector.FilterRad) * Perp ).sum().__abs__() * Detector.Meshes.dOmega.Radian )**2

    ParaFiltre = ( ( np.sin(Detector.FilterRad) * Para ).sum().__abs__() * Detector.Meshes.dOmega.Radian )**2

    return np.asscalar(PerpFiltre + ParaFiltre)


def Coupling_NoFilter(Detector, Scatterer):

    Perp = Detector.Fourier *\
           (Scatterer.Field.Perpendicular).__abs__() *\
           (np.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()


    Para = Detector.Fourier *\
           (Scatterer.Field.Parallel).__abs__() *\
           (np.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()

    PerpFiltre = (Perp.sum().__abs__() * Detector.Meshes.dOmega.Radian )**2

    ParaFiltre = (Para.sum().__abs__() * Detector.Meshes.dOmega.Radian )**2

    return np.asscalar(PerpFiltre + ParaFiltre)



def Coupling(Scatterer, Detector, Polarization):

    if Polarization in ['Parallel']:
        return Coupling_Para(Detector, Scatterer)


    if Polarization in ['Perpendicular']:
        return Coupling_Perp(Detector, Scatterer)


    if Polarization in ['Filtered']:

        return Coupling_Filtered(Detector, Scatterer)


    if Polarization in ['NoFiltered']:

        return Coupling_NoFilter(Detector, Scatterer)


    if Polarization in ['all']:
        CouplingDict = {}

        coupling = Coupling_Para(Detector, Scatterer)
        CouplingDict['Parallel'] = coupling

        coupling = Coupling_Perp(Detector, Scatterer)
        CouplingDict['Perpendicular'] = coupling

        coupling = Coupling_Filtered(Detector, Scatterer)
        CouplingDict['Filtered'] = coupling

        coupling = Coupling_NoFilter(Detector, Scatterer)
        CouplingDict['NoFiltered'] = coupling

        return CouplingDict
