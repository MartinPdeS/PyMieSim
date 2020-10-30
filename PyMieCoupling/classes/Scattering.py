from numpy import pi
import numpy as np

import scipy.interpolate as interp
import matplotlib.pyplot as plt
from matplotlib import cm
import PyMieScatt
import functools

from PyMieCoupling.classes.Fields import Field
from PyMieCoupling.classes.Meshes import Meshes as MieMesh
from PyMieCoupling.classes.Plots import S1S2Plot
from PyMieCoupling.classes.Representations import S1S2
from PyMieCoupling.classes.Misc import Source, Make3D, S1S2ToField
from typing import Tuple

global i
i = complex(0, 1)


class Scatterer(object):
    """Object containing all scatterer-related attributes.

    Parameters
    ----------
    diameter : float
        Diameter of the scatterer.
    wavelength : float
        Wavelength of the incident lightfield.
    index : float
        Refractive index of the scatterer.
    npts : int
        Number of points for the full solid angle of the far-field, later to
        be interpolated.

    Attributes
    ----------
    Full : <Fields class>
        It represents the entire Far-field representation of the scatterer.
    ComputeS1S2 : type
        Methode using package PyMieScatt to compute S1 and S2 parameter form mu value.
    diameter
    wavelength
    index
    npts

    """

    def __init__(self,
                 Diameter: float,
                 Source: Source,
                 Index: float,
                 Npts: int = None,
                 Meshes: MieMesh = None,
                 ThetaBound: list = [-180, 180],
                 ThetaOffset: float = 0,
                 PhiBound: list = [-180, 180],
                 PhiOffset: float = 0,
                 CacheTrunk: int = 0) -> None:

        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter * 2 )

        self.CacheTrunk = CacheTrunk

        self.GenMesh(Meshes, ThetaBound, PhiBound, ThetaOffset, PhiOffset, Npts)

        self.__S1S2 = None

        self.GenField()


    def GenMesh(self,
                Meshes: MieMesh = None,
                ThetaBound: list = [-90, 90],
                PhiBound: list = [-90, 90],
                ThetaOffset: float = 0,
                PhiOffset: float = 0,
                Npts: int = 101):

        if Meshes:
            self.Meshes = Meshes
            assert not all([ThetaBound, PhiBound, ThetaOffset, PhiOffset, Npts])

        else:
            self.Meshes = MieMesh(ThetaBound = np.array(ThetaBound) + ThetaOffset,
                                  PhiBound   = np.array(PhiBound) + PhiOffset,
                                  Npts       = Npts)

        self.MuList = np.cos(self.Meshes.Phi.Vector.Radian)


    @property
    def S1S2(self) -> None:
        if self.__S1S2 is None:
            self.__S1S2 = S1S2(self.SizeParam,
                               self.MuList,
                               self.Index,
                               self.Meshes)
            return self.__S1S2

        else:
            return self.__S1S2



    def GenField(self):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.

        """

        Parallel, Perpendicular = S1S2ToField(S1S2         = self.S1S2,
                                              Source       = self.Source,
                                              Meshes       = self.Meshes)

        self.Field = Field(Perpendicular  = Perpendicular,
                           Parallel       = Parallel,
                           Meshes         = self.Meshes
                           )

    @property
    def Stokes(self) -> None:
        return self.Field.Stokes

    @property
    def Jones(self) -> None:
        return self.Field.Jones

    @property
    def SPF(self) -> None:
        return self.Field.SPF



    # -
