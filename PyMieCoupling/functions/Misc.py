import numpy as np
from typing import Tuple
import functools
import PyMieScatt
import cupy as cp
import fibermodes
from PyMieCoupling.classes.Misc import Operation as Op
from PyMieCoupling.functions.MieComputing import MieS1S2


def Make3D(item:      np.array,
           PhiMesh:   np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * Op.sin(PhiMesh) * Op.cos(ThetaMesh)

    Y = item * Op.sin(PhiMesh) * Op.sin(ThetaMesh)

    Z = item * Op.cos(PhiMesh)

    return X, Y, Z


def GetSPF(cuda=False, **kwargs):
    return ComputeSPF(**kwargs)



def GetStokes(Parallel:      cp.ndarray,
              Perpendicular: cp.ndarray,
              cuda:          bool       = False) -> cp.ndarray:

    Array = Op.empty(cuda)( [4, *Parallel.shape] )

    I = Parallel.__abs__()**2 + Perpendicular.__abs__()**2

    Array[0,:,:] = I

    Array[1,:,:] = (Parallel.__abs__()**2 - Perpendicular.__abs__()**2)/I

    Array[2,:,:] = 2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    Array[3,:,:] = -2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    return Array


def GetJones(Parallel:      cp.ndarray,
             Perpendicular: cp.ndarray,
             cuda:          bool       = False) -> cp.ndarray:

    Array = Op.empty(cuda)( [2, * Parallel.shape] )

    delta = Op.angle(cuda)(Parallel) - Op.angle(cuda)(Perpendicular)

    A = Parallel.__abs__() / Op.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    B = Perpendicular.__abs__() / Op.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    return Op.array(cuda)([A, B * Op.exp(complex(0,1)*delta)])



def ComputeSPF(Parallel: np.ndarray, Perpendicular: np.ndarray) -> np.ndarray:
    return Parallel.__abs__()**2 + Perpendicular.__abs__()**2



def GetS1S2(Index,
            SizeParam,
            Meshes,
            cuda,
            CacheTrunk=None) -> Tuple[list, list]:

    MuList = np.cos(Meshes.Phi.Vector.Radian.tolist())

    if CacheTrunk: MuList = Op.round(cuda)(MuList, CacheTrunk)

    S1, S2 = [], []

    for Mu in MuList:

        temp0, temp1 = MieS1S2(Index,
                                SizeParam,
                                Mu)

        S1.append(temp0)
        S2.append(temp1)

    return Op.array(cuda)([S1, S2])



def _GetS1S2(Index,
                SizeParam,
                Meshes,
                cuda,
                CacheTrunk=None) -> Tuple[list, list]:

    MuList = np.cos(Meshes.Phi.Vector.Radian.tolist())

    if CacheTrunk: MuList = Op.round(cuda)(MuList, CacheTrunk)

    S1, S2 = np.ones(len(MuList)), np.ones(len(MuList))

    for nm, Mu in enumerate(MuList):

        S1[nm], S2[nm] = PyMieScatt.MieS1S2(Index,
                                    SizeParam,
                                    Mu)


    return Op.array(cuda)([S1, S2])





@functools.lru_cache(maxsize=201)
def WrapS1S2(Mu, Index, SizeParam) -> Tuple[float, float]:

    S1, S2 = PyMieScatt.MieS1S2(Index,
                                SizeParam,
                                Mu)

    return S1, S2





def S1S2ToField(S1S2, Source, Meshes, cuda):

    if Source.Polarization is not None:

        Parallel = Op.outer(S1S2.Array[0], Op.sin(Meshes.Theta.Vector.Radian))

        Perpendicular = Op.outer(S1S2.Array[1], Op.cos(Meshes.Theta.Vector.Radian))

    else:

        Parallel = Op.outer( S1S2.Array[0],  Op.ones(cuda)( len( S1S2.Array[0] ) )/Op.sqrt(2) )

        Perpendicular = Op.outer( S1S2.Array[1], Op.ones(cuda)( ( S1S2.Array[1] ) )/Op.sqrt(2) )

    return Parallel, Perpendicular




def GetLP(Fiber,
          Mode,
          Wavelength: float,
          Size:       float,
          Npts:       int,
          cuda:       bool):

    Field = fibermodes.field.Field(Fiber,
                                   Mode,
                                   Wavelength,
                                   Size,
                                   Npts).Ex()

    Field = Op.array(cuda)(Field)

    Field /= (Field.__abs__()).sum()

    Fourier = Op.fft(cuda).fft2(Field)

    Fourier /= GenShift(Npts, cuda)

    Fourier = Op.fft(cuda).fftshift(Fourier)

    Fourier /= (Fourier.__abs__()).sum()

    return Field, Fourier


def GenShift(Npts, cuda):

    phase_shift = Op.exp(-complex(0, 1) * Op.pi * Op.arange(cuda)(Npts)*(Npts-1)/Npts)

    shift_grid, _ = Op.meshgrid(cuda)(phase_shift, phase_shift)

    return shift_grid * shift_grid.T




# -
