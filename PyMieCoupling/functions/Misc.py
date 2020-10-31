import numpy as np
from typing import Tuple
import functools
import PyMieScatt


def Make3D(item: np.array,
           PhiMesh: np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z




def ComputeStokes(Parallel, Perpendicular):

    Array = np.empty( [4, *np.shape(Parallel)] )

    I = np.abs(Parallel)**2 + np.abs(Perpendicular)**2

    Array[0,:,:] = I

    Array[1,:,:] = (np.abs(Parallel)**2 - np.abs(Perpendicular)**2)/I

    Array[2,:,:] = 2*np.real( Parallel * np.conjugate(Perpendicular))/I

    Array[3,:,:] = -2*np.imag( Parallel * np.conjugate(Perpendicular))/I

    return Array



def ComputeJones(Parallel, Perpendicular):

    Array = np.empty( [2, *np.shape(Parallel)] )

    delta = np.angle(Parallel)-np.angle(Perpendicular)

    A = np.abs(Parallel) / np.sqrt(abs(Parallel)**2 + np.abs(Perpendicular)**2)

    B = np.abs(Perpendicular) / np.sqrt(abs(Parallel)**2 + np.abs(Perpendicular)**2)

    return np.array([A, B*np.exp(complex(0,1)*delta)])



def ComputeSPF(Parallel, Perpendicular):

    return np.abs(Parallel)**2 + np.abs(Perpendicular)**2





def ComputeS1S2(MuList, Index, SizeParam, CacheTrunk=None):

    if CacheTrunk: MuList = np.round(MuList, CacheTrunk)

    S1, S2 = [], []

    for Mu in MuList:

        temp0, temp1 = WrapS1S2(Mu, Index, SizeParam)

        S1.append(temp0)
        S2.append(temp1)

    return [S1, S2]


@functools.lru_cache(maxsize=201)
def WrapS1S2(Mu, Index, SizeParam) -> Tuple[float, float]:

    S1, S2 = PyMieScatt.MieS1S2(Index,
                                SizeParam,
                                Mu)

    return S1, S2















# -
