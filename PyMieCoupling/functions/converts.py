import numpy as np
from typing import Union

try:
    import cupy as cp
except:
    import numpy as cp


def rad2deg(RadSpace) -> Union[cp.ndarray, np.ndarray]:

        return RadSpace * 180/np.pi


def deg2rad(AngleSpace) -> Union[cp.ndarray, np.ndarray]:

    return AngleSpace * np.pi/180



def Angle2Direct(AngleVec: np.ndarray,
                 k:        float,
                 cuda:     bool = False) -> Union[cp.ndarray, np.ndarray]:

    RadSpace = AngleVec * np.pi / 180

    FourierSpace = np.sin(RadSpace) * k / (2 * np.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = np.fft.fftshift( np.fft.fftfreq( AngleVec.shape[0], d = fourier_unit ) )

    return DirectSpace


def Direct2Angle(DirectVec: np.ndarray,
                 k:         float,
                 cuda:      bool = False) -> Union[cp.ndarray, np.ndarray]:

    direct_unit = (DirectVec[1] - DirectVec[0]).__abs__()

    FourierSpace = np.fft.fftshift( np.fft.fftfreq( DirectVec.shape[0], d = direct_unit ) )

    AngleVec = np.arcsin(2 * np.pi * FourierSpace / k) #conversion spatial frequency to angular space

    if np.nan in AngleVec or cp.nan in AngleVec:
        raise Exception('Warning ill defined resolution -> angle definition!')

    return AngleVec * 180 / np.pi



def NA2Angle(NA:   float,
             cuda: bool = False) -> Union[cp.ndarray, np.ndarray]:

    Angle = rad2deg( np.arcsin(NA) )

    __ThetaBound = np.array( [-Angle, Angle] )

    __PhiBound = np.array( [-Angle, Angle] )


    return __ThetaBound, __PhiBound



def CuPy2NumPy(*items):
    ItemList = []
    for item in items:
        if isinstance(item, np.ndarray):
            ItemList.append(item)
        else:
            ItemList.append( cp.asnumpy(item) )

    temp = tuple(ItemList)

    return (*temp,)






# -
