import numpy as np
from ai import cs

from PyMieSim.Physics import Angle

""" Ref: https://optiwave.com/optifdtd-manuals/fdtd-far-field-transform/"""

def rad2deg(RadSpace) -> np.ndarray:

        return RadSpace * 180/np.pi



def deg2rad(AngleSpace) -> np.ndarray:

    return AngleSpace * np.pi/180



def Angle2Direct(AngleVec: np.ndarray, k: float,) -> np.ndarray:

    RadSpace = np.deg2rad(AngleVec)

    FourierSpace = np.sin(RadSpace) * k / (2 * np.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = np.fft.fftshift( np.fft.fftfreq( AngleVec.shape[0], d = fourier_unit ) )

    return DirectSpace



def Direct2Angle(DirectVec: np.ndarray, k: float) -> np.ndarray:

    direct_unit = (DirectVec[1] - DirectVec[0]).__abs__()

    FourierSpace = np.fft.fftshift( np.fft.fftfreq( DirectVec.shape[0], d = direct_unit ) )

    AngleVec = np.arcsin(2 * np.pi * FourierSpace / k) # conversion spatial frequency to angular space

    if np.isnan(AngleVec).any():
        raise Exception("Magnification too large.")

    return AngleVec * 180 / np.pi



def NA2Angle(NA: float) -> np.ndarray:

    return Angle( np.arcsin(NA), unit='Radian')



def Direct2spherical(X, Y, MaxAngle):
    Z = 50 / np.tan(MaxAngle)

    _, Phi, Theta = cs.cart2sp(X, Y, X*0+Z)

    return Phi, Theta

def Direct2Angle(X, Y, MaxAngle):
    MaxZ = np.max(X) / np.cos(MaxAngle)



def AngleUnit2DirectUnit(Angle, k):
    FourierSpace = np.sin(Angle) * k / (2 * np.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = np.fft.fftshift( np.fft.fftfreq( Angle.shape[0], d = fourier_unit ) )

    return DirectSpace
# -
