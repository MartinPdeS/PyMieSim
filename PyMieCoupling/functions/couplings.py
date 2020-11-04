import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Meshes import Meshes as MieMesh
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Misc import Operation as Op
from typing import Union

def PointFieldCoupling(Detector: Union[LPmode, Photodiode],
                       Source: Scatterer):

    dOmega = Detector.Meshes.Phi.Delta.Radian *\
             Detector.Meshes.Theta.Delta.Radian

    if Detector._coupling == 'Amplitude':

        Perp = Detector.Field.Array *\
               Source.Field.Perpendicular *\
               (Op.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()

        Perp = ( Perp.sum() ).__abs__()**2

        Para = Detector.Field.Array *\
               Source.Field.Parallel *\
               (Op.sin(Detector.Meshes.Phi.Mesh.Radian).T).__abs__()

        Para = (Para * dOmega).sum().__abs__()**2


    elif Detector._coupling == 'Intensity':
        Perp = Detector.Fourier.Array *\
               (Source.Field.Perpendicular).__abs__() *\
               (Op.sin(Detector.Meshes.Phi.Mesh.Radian + Op.pi/2).T).__abs__()

        Perp = (Perp * dOmega).sum()**2

        Para = Detector.Fourier.Array *\
               (Source.Field.Parallel).__abs__() *\
               (Op.sin(Detector.Meshes.Phi.Mesh.Radian + Op.pi/2).T).__abs__()

        Para = (Para * dOmega).sum()**2

    return Para, Perp




def MeanFieldCoupling(Field0: np.array,
                      Field1: np.array):

    temp = Field0 * Field1

    temp = np.fft.fftshift(temp)

    temp = np.fft.ifft2(temp)

    temp = np.fft.fftshift(temp)

    return (temp).__abs__()**2
