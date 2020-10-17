from numpy import pi
import numpy as np

import scipy.interpolate as interp
import matplotlib.pyplot as plt
from matplotlib import cm
import PyMieScatt

from miecoupling.src.functions.converts import rad2deg, deg2rad
from miecoupling.src.classes.Fields import Field
from miecoupling.src.classes.Meshes import Meshes
from miecoupling.src.classes.Plots import S1S2Plot

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
                 npts: int = 201,
                 ThetaBound: list = [-180, 180],
                 PhiBound: list = [-180, 180]):

        self.diameter = diameter

        self.wavelength = wavelength

        self.index = index

        self.npts = npts

        self.Full = Meshes(ThetaBound=ThetaBound,
                           PhiBound=PhiBound,
                           npts=npts)

        self.computeS1S2()


    def PlotFields(self):

        SPF = self.Make3DField()

        Plot = S1S2Plot(np.abs(self.S1), np.abs(self.S2), *SPF, self.Full)


    def Make3DField(self,):

        Theta = np.linspace(0, 2*np.pi, self.npts)

        Phi = np.linspace(0, np.pi, self.npts)

        FuncSPF = interp.interp2d(self.Field.ThetaVec.Radian,
                                  self.Field.PhiVec.Radian,
                                  self.Field.SPF,
                                  kind='cubic')


        SPF = FuncSPF(self.Field.ThetaVec.Radian, self.Field.PhiVec.Radian)


        X = SPF * np.sin(self.Field.PhiMesh.Radian) * np.cos(self.Field.ThetaMesh.Radian)
        Y = SPF * np.sin(self.Field.PhiMesh.Radian) * np.sin(self.Field.ThetaMesh.Radian)
        Z = SPF * np.cos(self.Field.PhiMesh.Radian)

        return [X, Y, Z]


    def computeS1S2(self):

        MuList = np.cos(self.Full.PhiVec.Radian)

        self.S1, self.S2 = [], []

        for mu in MuList:

            SizeParam = np.pi*self.diameter/self.wavelength

            S1, S2 = PyMieScatt.MieS1S2(self.index,
                                        SizeParam,
                                        mu)

            self.S1.append(S1)
            self.S2.append(S2)



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

            Parallel = np.outer(self.S1, np.sin(self.Full.ThetaVec.Radian+np.pi/2))

            Perpendicular = np.outer(self.S2, np.cos(self.Full.ThetaVec.Radian+np.pi/2))

        elif PolarizationAngle is None:

            self.Polarization = "None"

            Parallel = np.outer(self.S1,  np.ones(self.npts)/np.sqrt(2))

            Perpendicular = np.outer(self.S2, np.ones(self.npts)/np.sqrt(2))


        self.Field = Field(Perpendicular=Perpendicular,
                           Parallel=Parallel,
                           Meshes=self.Full
                           )



    # -
