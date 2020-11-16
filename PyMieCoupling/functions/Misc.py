import numpy as np
from typing import Tuple
import functools
import fibermodes
from PyMieCoupling.cpp.S1S2 import MieS1S2



def Make3D(item:      np.array,
           PhiMesh:   np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z


def GetSPF(**kwargs):
    return ComputeSPF(**kwargs)



def GetStokes(Parallel:      np.ndarray,
              Perpendicular: np.ndarray) -> np.ndarray:

    Array = np.empty( [4, *Parallel.shape] )

    I = Parallel.__abs__()**2 + Perpendicular.__abs__()**2

    Array[0,:,:] = I

    Array[1,:,:] = (Parallel.__abs__()**2 - Perpendicular.__abs__()**2)/I

    Array[2,:,:] = 2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    Array[3,:,:] = -2 * ( Parallel * Perpendicular.conjugate() ).imag / I

    return Array


def GetJones(Parallel:      np.ndarray,
             Perpendicular: np.ndarray) -> np.ndarray:

    Array = np.empty( [2, * Parallel.shape] )

    delta = np.angle(Parallel) - np.angle(Perpendicular)

    A = Parallel.__abs__() / np.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    B = Perpendicular.__abs__() / np.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    return np.array([A, B * np.exp(complex(0,1)*delta)])



def ComputeSPF(Parallel: np.ndarray, Perpendicular: np.ndarray) -> np.ndarray:
    return Parallel.__abs__()**2 + Perpendicular.__abs__()**2



def GetS1S2(Index,
            SizeParam,
            Meshes) -> Tuple[np.ndarray, np.ndarray]:

    S1, S2 = MieS1S2(Index,
                     SizeParam,
                     Meshes.Phi.Vector.Radian.tolist(),
                     Meshes.Theta.Vector.Radian.tolist(),
                     )

    return np.array(S1), np.array(S2)





def S1S2ToField(S1,
                S2,
                Meshes) -> Tuple[np.ndarray, np.ndarray]:


    Parallel = np.outer(S1, np.sin(Meshes.Theta.Vector.Radian))

    Perpendicular = np.outer(S2, np.cos(Meshes.Theta.Vector.Radian))


    return Parallel, Perpendicular





def GetLP(Fiber,
          Mode,
          Wavelength: float,
          Size:       float,
          Npts:       int):

    Field = fibermodes.field.Field(Fiber,
                                   Mode,
                                   Wavelength,
                                   Size,
                                   Npts).Ex()

    Field = np.array(Field)

    Field /= (Field.__abs__()).sum()

    Fourier = np.fft.fft2(Field)

    Fourier /= GenShift(Npts)

    Fourier = np.fft.fftshift(Fourier)

    Fourier /= (Fourier.__abs__()).sum()

    return Field, Fourier


def GenShift(Npts):

    phase_shift = np.exp(-complex(0, 1) * np.pi * np.arange(Npts)*(Npts-1)/Npts)

    shift_grid, _ = np.meshgrid(phase_shift, phase_shift)

    return shift_grid * shift_grid.T



def Coupling(Field:         np.ndarray,
             Parallel:      np.ndarray,
             Perpendicular: np.ndarray,
             Phi:           np.ndarray):


    dOmega = (Phi[0,0] - Phi[0,1]).__abs__()**2

    Perp = Field * Perpendicular * np.sin(Phi).__abs__()

    Perp = ( Perp * dOmega).sum().__abs__()**2

    Para = Field * Parallel * np.sin(Phi).__abs__()

    Para = (Para * dOmega).sum().__abs__()**2

    return Para, Perp


# -
