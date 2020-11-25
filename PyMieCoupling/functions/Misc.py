import numpy as np
from typing import Tuple
import functools
import fibermodes
from PyMieCoupling.cpp.S1S2 import MieS1S2







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




# -
