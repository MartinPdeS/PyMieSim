from numpy import pi
import numpy as np

import scipy.interpolate as interp
import matplotlib.pyplot as plt
from matplotlib import cm
import PyMieScatt
import functools

from PyMieCoupling.functions.converts import rad2deg, deg2rad
from PyMieCoupling.classes.Fields import Field
from PyMieCoupling.classes.Meshes import Meshes as MieMesh
from PyMieCoupling.classes.Plots import S1S2Plot

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
    computeS1S2 : type
        Methode using package PyMieScatt to compute S1 and S2 parameter form mu value.
    diameter
    wavelength
    index
    npts

    """

    def __init__(self,
                 diameter: float,
                 wavelength: float,
                 index: float,
                 npts: int = None,
                 Meshes: MieMesh = None,
                 ThetaBound: list = [-180, 180],
                 ThetaOffset: float = 0,
                 PhiBound: list = [-180, 180],
                 PhiOffset: float = 0,
                 CacheTrunk: int = 0):

        self.diameter, self.wavelength = diameter, wavelength

        self.index= index

        self.CacheTrunk = CacheTrunk

        if Meshes:
            self.Meshes = Meshes
            assert not all([ThetaBound, PhiBound, ThetaOffset, PhiOffset, npts])

        else:
            self.Meshes = MieMesh(ThetaBound = np.array(ThetaBound) + ThetaOffset,
                                  PhiBound   = np.array(PhiBound) + PhiOffset,
                                  npts       = npts)

        self.computeS1S2()

        self.GenField(PolarizationAngle=0)


    def PlotS1S2(self):

        SPF = self.Make3DField(self.Field.SPF)

        Plot = S1S2Plot(np.abs(self.S1), np.abs(self.S2), *SPF, self.Meshes)


    def Make3DField(self, item):

        X = item * np.sin(self.Field.Phi.Mesh.Radian) * np.cos(self.Field.Theta.Mesh.Radian)

        Y = item * np.sin(self.Field.Phi.Mesh.Radian) * np.sin(self.Field.Theta.Mesh.Radian)

        Z = item * np.cos(self.Field.Phi.Mesh.Radian)

        return [X, Y, Z]


    def computeS1S2(self):

        MuList = np.cos(self.Meshes.Phi.Vector.Radian)

        if self.CacheTrunk: MuList = np.round(MuList, self.CacheTrunk)

        self.S1, self.S2 = [], []

        for Mu in MuList:

            self.SizeParam = np.pi * self.diameter / self.wavelength

            S1, S2 = self.WrapS1S2(Mu)

            self.S1.append(S1)
            self.S2.append(S2)


    @functools.lru_cache(maxsize=201)
    def WrapS1S2(self, Mu):

        S1, S2 = PyMieScatt.MieS1S2(self.index,
                                    self.SizeParam,
                                    Mu)

        return S1, S2


    def GenField(self, PolarizationAngle=0):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.

        Parameters
        ----------
        PolarizationAngle : float
            Value representing the polarization angle of the inciendent field.

        """

        if PolarizationAngle is not None:
            Polarization = deg2rad(PolarizationAngle)

            self.Polarization = Polarization

            Parallel = np.outer(self.S1, np.sin(self.Meshes.Theta.Vector.Radian))

            Perpendicular = np.outer(self.S2, np.cos(self.Meshes.Theta.Vector.Radian))

        elif PolarizationAngle is None:

            self.Polarization = "None"

            Parallel = np.outer(self.S1,  np.ones(len(self.S1))/np.sqrt(2))

            Perpendicular = np.outer(self.S2, np.ones((self.S2))/np.sqrt(2))


        self.Field = Field(Perpendicular  = Perpendicular,
                           Parallel       = Parallel,
                           Meshes         = self.Meshes
                           )



    # -
