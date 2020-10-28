import numpy as np
from typing import Tuple


def Make3D(item: np.array,
           PhiMesh: np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z
