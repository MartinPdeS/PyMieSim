import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import Detector
from PyMieCoupling.classes.Meshes import Meshes as MieMesh
from PyMieCoupling.classes.Scattering import Scatterer


def PointFieldCoupling(Detector: Detector,
                       Source: Scatterer,
                       Field: str = None):

    if not Field:
        raise Exception('Field must be specified [Parallel, Perpendicular]')

    if Field == 'Parallel':
        Source = Source.Field.Parallel

    elif Field == 'Perpendicular':
        Source = Source.Field.Perpendicular

    if Detector._coupling == 'Amplitude':

        temp = Detector.Field * Source * np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian).T)

        temp = np.sum(temp)

        temp = np.abs(temp)**2

        return temp


    elif Detector._coupling == 'Intensity':

        temp = np.abs(Source) * Detector.Field * np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian).T)

        temp = np.sum(temp)

        return temp




def MeanFieldCoupling(Field0: np.array,
                      Field1: np.array):

    temp = Field0 * Field1

    temp = np.fft.fftshift(temp)

    temp = np.fft.ifft2(temp)

    temp = np.fft.fftshift(temp)

    return np.abs(temp)**2
