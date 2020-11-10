import numpy as np
from PyMieCoupling.classes.Misc import Operation as Op
from typing import Union

try:
    import cupy as cp
except:
    import numpy as cp


def rad2deg(RadSpace) -> Union[cp.ndarray, np.ndarray]:

        return RadSpace * 180/Op.pi


def deg2rad(AngleSpace) -> Union[cp.ndarray, np.ndarray]:

    return AngleSpace * Op.pi/180



def Angle2Direct(AngleVec: np.ndarray,
                 k:        float,
                 cuda:     bool = False) -> Union[cp.ndarray, np.ndarray]:

    RadSpace = AngleVec * Op.pi / 180

    FourierSpace = Op.sin(RadSpace) * k / (2 * Op.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = Op.fft(cuda).fftshift( Op.fft.fftfreq( AngleVec.shape[0], d = fourier_unit ) )

    return DirectSpace


def Direct2Angle(DirectVec: np.ndarray,
                 k:         float,
                 cuda:      bool = False) -> Union[cp.ndarray, np.ndarray]:

    direct_unit = (DirectVec[1] - DirectVec[0]).__abs__()

    FourierSpace = Op.fft(cuda).fftshift( Op.fft(cuda).fftfreq( DirectVec.shape[0], d = direct_unit ) )

    AngleVec = Op.arcsin(2 * Op.pi * FourierSpace / k) #conversion spatial frequency to angular space

    if np.nan in AngleVec or cp.nan in AngleVec:
        raise Exception('Warning ill defined resolution -> angle definition!')

    return AngleVec * 180 / Op.pi



def NA2Angle(NA:   float,
             cuda: bool = False) -> Union[cp.ndarray, np.ndarray]:

    Angle = rad2deg( Op.arcsin(NA) )

    __ThetaBound = Op.array(cuda)( [-Angle, Angle] )

    __PhiBound = Op.array(cuda)( [-Angle, Angle] )


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
