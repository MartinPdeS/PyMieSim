import numpy as np
import matplotlib.pyplot as plt


def PointFieldCoupling(Detector,
                       Source,
                       Mesh):


    if Detector._coupling == 'Amplitude':

        temp = Detector.Field * Source * np.sin(Mesh.PhiMesh.Radian+np.pi/2)

        temp = np.sum(temp)**2

        temp = np.abs(temp)

        return temp

    elif Detector._coupling == 'Intensity':

        temp = np.abs(Source) * Detector.Field * np.sin(Mesh.PhiMesh.Radian)

        temp = np.sum(temp)

        return temp




def MeanFieldCoupling(Field0, Field1):

    temp = Field0 * Field1

    temp = np.fft.fftshift(temp)

    temp = np.fft.ifft2(temp)

    temp = np.fft.fftshift(temp)

    return np.abs(temp)**2


def PointIntensityCoupling():
    pass


def MeanIntensityCoupling():
    pass
