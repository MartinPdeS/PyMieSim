import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import DetectorMeta
from PyMieCoupling.classes.Meshes import Meshes as MieMesh
from PyMieCoupling.classes.Scattering import Scatterer


def PointFieldCoupling(Detector: DetectorMeta,
                       Source: Scatterer):

    dOmega = Detector.Meshes.Phi.Delta.Radian *\
             Detector.Meshes.Theta.Delta.Radian

    if Detector._coupling == 'Amplitude':

        Perp = Detector.Field *\
               Source.Field.Perpendicular *\
               np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian).T)

        Perp = np.abs( np.sum(Perp) )**2

        Para = Detector.Field *\
               Source.Field.Parallel *\
               np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian).T)

        Para = np.abs( np.sum(Para * dOmega) )**2


    elif Detector._coupling == 'Intensity':

        Perp = Detector.Fourier *\
               np.abs(Source.Field.Perpendicular) *\
               np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian).T)

        Perp = np.sum(Perp * dOmega)**2

        Para = Detector.Fourier *\
               np.abs(Source.Field.Parallel) *\
               np.abs(np.sin(Detector.Meshes.Phi.Mesh.Radian).T)

        Para = np.sum(Para * dOmega)**2

        data = np.sum(Source.Field.Parallel) / dOmega


    return Para, Perp




def MeanFieldCoupling(Field0: np.array,
                      Field1: np.array):

    temp = Field0 * Field1

    temp = np.fft.fftshift(temp)

    temp = np.fft.ifft2(temp)

    temp = np.fft.fftshift(temp)

    return np.abs(temp)**2
