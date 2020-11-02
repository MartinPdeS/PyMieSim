import numpy as np
import cupy as cp

def rad2deg(RadSpace):

    if isinstance(RadSpace, np.ndarray) or isinstance(RadSpace, np.float):
        return RadSpace * 180/np.pi

    if isinstance(RadSpace, cp.ndarray) or isinstance(RadSpace, cp.float):
        return RadSpace * 180/cp.pi


def deg2rad(AngleSpace):
    if isinstance(AngleSpace, np.ndarray) or isinstance(AngleSpace, np.float):
        return AngleSpace * np.pi/180

    if isinstance(AngleSpace, cp.ndarray) or isinstance(AngleSpace, cp.float):
        return AngleSpace * cp.pi/180


def Angle2Direct(AngleVec, k):

    if isinstance(RadSpace, np.ndarray) or isinstance(RadSpace, np.float):
        RadSpace = AngleVec*np.pi/180

        FourierSpace = np.sin(RadSpace)*k/(2*np.pi)

        fourier_unit = (FourierSpace[1]-FourierSpace[0]).__abs__()

        DirectSpace = np.fft.fftshift(np.fft.fftfreq(np.shape(AngleVec)[0], d=fourier_unit))

    if isinstance(RadSpace, cp.ndarray) or isinstance(RadSpace, cp.float):
        RadSpace = AngleVec*cp.pi/180

        FourierSpace = cp.sin(RadSpace)*k/(2*cp.pi)

        fourier_unit = (FourierSpace[1]-FourierSpace[0]).__abs__()

        DirectSpace = cp.fft.fftshift(cp.fft.fftfreq(cp.shape(AngleVec)[0], d=fourier_unit))

    return DirectSpace


def Direct2Angle(DirectVec, k):

    if isinstance(DirectVec, np.ndarray) or isinstance(DirectVec, np.float):

        direct_unit = (DirectVec[1]-DirectVec[0]).__abs__()

        FourierSpace = np.fft.fftshift(np.fft.fftfreq(np.shape(DirectVec)[0],d=direct_unit))

        AngleVec = np.arcsin(2*np.pi*FourierSpace/ k) #conversion spatial frequency to angular space

        if np.nan in AngleVec:
            raise Exception('Warning ill defined resolution -> angle definition!')

        AngleVec = AngleVec* 180/np.pi


    if isinstance(DirectVec, cp.ndarray) or isinstance(DirectVec, cp.float):

        direct_unit = (DirectVec[1]-DirectVec[0]).__abs__()

        FourierSpace = cp.fft.fftshift(cp.fft.fftfreq(cp.shape(DirectVec)[0],d=direct_unit))

        AngleVec = cp.arcsin(2*cp.pi*FourierSpace/ k) #conversion spatial frequency to angular space

        AngleVec = AngleVec* 180/cp.pi

        if cp.nan in AngleVec:
            raise Exception('Error in angle definition!')

    return AngleVec



def NA2Angle(NA, GPU):

    if GPU:

        Angle = rad2deg( cp.arcsin(NA) )

        __ThetaBound = cp.array([-Angle, Angle])

        __PhiBound = cp.array([-Angle, Angle])

    else:

        Angle = rad2deg( np.arcsin(NA) )

        __ThetaBound = np.array([-Angle, Angle])

        __PhiBound = np.array([-Angle, Angle])

    return __ThetaBound, __PhiBound



def CuPy2NumPy(*items):
    ItemList = []
    for item in items:
        if isinstance(item, cp.ndarray):
            ItemList.append(cp.asnumpy(item))
        else:
            ItemList.append(item)

    temp = tuple(ItemList)

    return (*temp,)






# -
