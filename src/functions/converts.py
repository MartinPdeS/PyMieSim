import numpy as np


def rad2deg(RadSpace):
    return RadSpace * 180/np.pi


def deg2rad(AngleSpace):
    return AngleSpace * np.pi/180

def Angle2Direct(AngleVec, k):

    RadSpace = AngleVec*np.pi/180

    FourierSpace = np.sin(RadSpace)*k/(2*np.pi)

    fourier_unit = np.abs(FourierSpace[1]-FourierSpace[0])

    DirectSpace = np.fft.fftshift(np.fft.fftfreq(np.shape(AngleVec)[0], d=fourier_unit))

    return DirectSpace


def Direct2Angle(DirectVec, k):

    direct_unit = np.abs(DirectVec[1]-DirectVec[0])

    FourierSpace = np.fft.fftshift(np.fft.fftfreq(np.shape(DirectVec)[0],d=direct_unit))

    AngleVec = np.arcsin(2*np.pi*FourierSpace/ k) #conversion spatial frequency to angular space

    AngleVec = AngleVec* 180/np.pi

    return AngleVec
