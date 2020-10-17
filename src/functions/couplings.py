import numpy as np
import matplotlib.pyplot as plt


def PointFieldCoupling(Field0, Field1):

    Field0, Field1 = np.array(Field0), np.array(Field1)

    return np.abs(np.sum(Field0*Field1)**2)


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
