import numpy as np
from typing import Tuple
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patches as mpatches


def Make3D(item: np.array,
           PhiMesh: np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z
