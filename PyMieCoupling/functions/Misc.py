import numpy as np
from typing import Tuple
import functools
import PyMieScatt
import cupy as cp
import fibermodes


def Make3D(item: np.array,
           PhiMesh: np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z


def GetStokes(GPU=False, **kwargs):
    if GPU:
        return ComputeStokesGPU(**kwargs)
    else:
        return ComputeStokes(**kwargs)


def GetJones(GPU=False, **kwargs):
    if GPU:
        return ComputeJonesGPU(**kwargs)
    else:
        return ComputeJones(**kwargs)


def GetSPF(GPU=False, **kwargs):
    return ComputeSPF(**kwargs)


def GetS1S2(GPU=False, **kwargs):
    if GPU:
        return ComputeS1S2GPU(**kwargs)
    else:
        return ComputeS1S2(**kwargs)



def ComputeStokesGPU(Parallel: cp.ndarray, Perpendicular: cp.ndarray) -> cp.ndarray:

    Array = cp.empty( [4, *Parallel.shape] )

    I = Parallel.__abs__()**2 + Perpendicular.__abs__()**2

    Array[0,:,:] = I

    Array[1,:,:] = (Parallel.__abs__()**2 - Perpendicular.__abs__()**2)/I

    Array[2,:,:] = 2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    Array[3,:,:] = -2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    return Array


def ComputeStokes(Parallel: np.ndarray, Perpendicular: np.ndarray) -> np.ndarray:

    Array = np.empty( [4, *Parallel.shape] )

    I = Parallel.__abs__()**2 + Perpendicular.__abs__()**2

    Array[0,:,:] = I

    Array[1,:,:] = (Parallel.__abs__()**2 - Perpendicular.__abs__()**2)/I

    Array[2,:,:] = 2 * ( Parallel * Perpendicular.conjugate() ).real / I

    Array[3,:,:] = -2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    return Array



def ComputeJones(Parallel: np.ndarray, Perpendicular: np.ndarray) -> np.ndarray:

    Array = np.empty( [2, *Parallel.shape] )

    delta = np.angle(Parallel)-np.angle(Perpendicular)

    A = Parallel._abs__() / np.sqrt(Parallel._abs__()**2 + Perpendicular._abs__()**2)

    B = Perpendicular._abs__() / np.sqrt(Parallel._abs__()**2 + Perpendicular._abs__()**2)

    return np.array([A, B*np.exp(complex(0,1)*delta)])


def ComputeJonesGPU(Parallel: cp.ndarray, Perpendicular: cp.ndarray) -> cp.ndarray:

    Array = cp.empty( [2, * Parallel.shape] )

    delta = cp.angle(Parallel)-cp.angle(Perpendicular)

    A = Parallel.__abs__() / cp.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    B = Perpendicular.__abs__() / cp.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    return cp.array([A, B*cp.exp(complex(0,1)*delta)])



def ComputeSPF(Parallel: np.ndarray, Perpendicular: np.ndarray) -> np.ndarray:

    return Parallel.__abs__()**2 + Perpendicular.__abs__()**2




def ComputeS1S2(Index, SizeParam, Meshes, CacheTrunk=None) -> Tuple[list, list]:

    MuList = np.cos(Meshes.Phi.Vector.Radian)

    if CacheTrunk: MuList = np.round(MuList, CacheTrunk)

    S1, S2 = [], []

    for Mu in MuList:

        temp0, temp1 = WrapS1S2(Mu, Index, SizeParam)

        S1.append(temp0)
        S2.append(temp1)

    return np.array([S1, S2])


@functools.lru_cache(maxsize=201)
def WrapS1S2(Mu, Index, SizeParam) -> Tuple[float, float]:

    S1, S2 = PyMieScatt.MieS1S2(Index,
                                SizeParam,
                                Mu)

    return S1, S2


def ComputeS1S2GPU(Index, SizeParam, Meshes, CacheTrunk=None) -> Tuple[list, list]:

    MuList = np.cos(cp.asnumpy(Meshes.Phi.Vector.Radian))

    if CacheTrunk: MuList = np.round(MuList, CacheTrunk)

    S1, S2 = [], []

    for Mu in MuList:

        temp0, temp1 = WrapS1S2(Mu, Index, SizeParam)
        S1.append(temp0)
        S2.append(temp1)

    return cp.array([S1, S2])



def S1S2ToField(GPU, **kwargs):

    if GPU:
        return _S1S2ToFieldGPU(**kwargs)

    else:
        return _S1S2ToField(**kwargs)



def _S1S2ToField(S1S2, Source, Meshes):

    if Source.Polarization is not None:

        Parallel = np.outer(S1S2.Array[0], np.sin(Meshes.Theta.Vector.Radian))

        Perpendicular = np.outer(S1S2.Array[1], np.cos(Meshes.Theta.Vector.Radian))

    else:

        Parallel = np.outer( S1S2.Array[0],  np.ones( len( S1S2.Array[0] ) )/np.sqrt(2) )

        Perpendicular = np.outer( S1S2.Array[1], np.ones( ( S1S2.Array[1] ) )/np.sqrt(2) )

    return Parallel, Perpendicular


def _S1S2ToFieldGPU(S1S2, Source, Meshes):

    if Source.Polarization is not None:

        Parallel = cp.outer(S1S2.Array[0], cp.sin(Meshes.Theta.Vector.Radian))

        Perpendicular = cp.outer(S1S2.Array[1], cp.cos(Meshes.Theta.Vector.Radian))

    else:

        Parallel = cp.outer( S1S2.Array[0],  cp.ones( len( S1S2.Array[0] ) )/cp.sqrt(2) )

        Perpendicular = cp.outer( S1S2.Array[1], cp.ones( ( S1S2.Array[1] ) )/cp.sqrt(2) )

    return Parallel, Perpendicular



def GetLP(GPU, **kwargs):

    if GPU: return GetLPGPU(**kwargs)

    else: return GetLPCPU(**kwargs)



def GetLPCPU(Fiber, Mode, Wavelength, Size, Npts):

    Field = fibermodes.field.Field(Fiber,
                                   Mode,
                                   Wavelength,
                                   Size,
                                   Npts).Ex()

    Field /= np.sum(Field.__abs__())

    Fourier = np.fft.fft2(Field)

    Fourier /= GenShift(Npts)

    Fourier = np.fft.fftshift(Fourier)

    Fourier /= np.sum(Fourier.__abs__())

    return Field, Fourier


def GetLPGPU(Fiber, Mode, Wavelength, Size, Npts):

    Field = fibermodes.field.Field(Fiber,
                                   Mode,
                                   Wavelength,
                                   Size,
                                   Npts).Ex()

    Field = cp.array(Field)

    Field /= cp.sum(Field.__abs__())

    Fourier = cp.fft.fft2(Field)

    Fourier /= cp.array( GenShift(Npts) )

    Fourier = cp.fft.fftshift(Fourier)

    Fourier /= cp.sum(Fourier.__abs__())

    return Field, Fourier


def GenShift(Npts):

    phase_shift = np.exp(-complex(0, 1)*np.pi*np.arange(Npts)*(Npts-1)/Npts)

    shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

    return shift_grid * shift_grid.T




# -
