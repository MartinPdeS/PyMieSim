import numpy as np
from typing import Tuple
from PyMieCoupling.functions.converts import deg2rad



class Source(object):

    def __init__(self,
                 Wavelength: float,
                 Polarization: float):

        self.Wavelength = Wavelength

        self.Polarization = Polarization

        self.k = 2 * np.pi / Wavelength



def Make3D(item: np.array,
           PhiMesh: np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z



def S1S2ToField(S1S2, Source, Meshes):

    if Source.Polarization is not None:

        Parallel = np.outer(S1S2.Array[0], np.sin(Meshes.Theta.Vector.Radian))

        Perpendicular = np.outer(S1S2.Array[1], np.cos(Meshes.Theta.Vector.Radian))

    else:

        Parallel = np.outer( S1S2.Array[0],  np.ones( len( S1S2.Array[0] ) )/np.sqrt(2) )

        Perpendicular = np.outer( S1S2.Array[1], np.ones( ( S1S2.Array[1] ) )/np.sqrt(2) )

    return Parallel, Perpendicular
