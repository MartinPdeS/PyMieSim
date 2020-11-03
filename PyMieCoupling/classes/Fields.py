
import numpy as np
import cupy as cp
from PyMieCoupling.classes.Representations import Stokes, Jones, SPF
from PyMieCoupling.classes.Meshes import Meshes as MieMesh


class Field(object):

    def __init__(self,
                 Perpendicular: np.ndarray,
                 Parallel:      np.ndarray,
                 Meshes:        MieMesh,
                 cuda:          bool):
        """
        Source -- https://www.physlab.org/wp-content/uploads/2016/07/Ch6-BYUOpticsBook_2013.pdf

        """
        self.__dict__ = Meshes.__dict__.copy()

        self.cuda = cuda

        self.Perpendicular, self.Parallel = Perpendicular, Parallel

        self.Meshes = Meshes

        self.__Jones, self.__Stokes, self.__SPF, self.__Delay, self.__Total = (None,)*5


    @property
    def Total(self) -> None:
        if self.__Total is None:
            self.__Total = self.ComputeTotal()
            return self.__Total

        else:
            return self.__Total


    @property
    def Delay(self) -> None:
        if self.__Delay is None:
            self.__Delay = self.ComputeDelay()
            return self.__Delay

        else:
            return self.__Delay


    @property
    def SPF(self) -> None:
        if self.__SPF is None:
            self.__SPF = SPF(Parallel      = self.Parallel,
                             Perpendicular = self.Perpendicular,
                             Meshes        = self.Meshes,
                             cuda          = self.cuda)
            return self.__SPF

        else:
            return self.__SPF


    @property
    def Stokes(self) -> None:
        if self.__Stokes is None:
            self.__Stokes = Stokes(Parallel      = self.Parallel,
                                   Perpendicular = self.Perpendicular,
                                   Meshes        = self.Meshes,
                                   cuda           = self.cuda)
            return self.__Stokes

        else:
            return self.__Stokes


    @property
    def Jones(self) -> None:
        if self.__Jones is None:
            self.__Jones = Jones(Parallel      = self.Parallel,
                                 Perpendicular = self.Perpendicular,
                                 Meshes        = self.Meshes,
                                 cuda           = self.cuda)
            return self.__Jones

        else:
            return self.__Jones


    def ComputeTotal(self) -> None:
        return sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)# * exp(complex(0,1)*self.delta)


    def ComputeDelay(self) -> None:
        return arctan(abs(self.Parallel)/abs(self.Perpendicular))


    def ComputeSPF(self) -> None:
        return abs(self.Parallel)**2 + abs(self.Perpendicular)**2















# -
