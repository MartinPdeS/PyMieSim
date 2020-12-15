import numpy as np

""" Ref: https://optiwave.com/optifdtd-manuals/fdtd-far-field-transform/"""

def rad2deg(RadSpace) -> np.ndarray:

        return RadSpace * 180/np.pi



def deg2rad(AngleSpace) -> np.ndarray:

    return AngleSpace * np.pi/180



def Angle2Direct(AngleVec: np.ndarray, k: float,) -> np.ndarray:

    RadSpace = AngleVec * np.pi / 180

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

    return rad2deg( np.arcsin(NA) ).squeeze()










# -
