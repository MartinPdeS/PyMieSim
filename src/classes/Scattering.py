from numpy import pi
import numpy as np

import scipy.interpolate as interp
import matplotlib.pyplot as plt
from matplotlib import cm
import PyMieScatt

from src.functions.converts import rad2deg, deg2rad
from src.classes.Fields import Field
from src.classes.Meshes import Meshes
from src.classes.Plots import S1S2Plot

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

    def __init__(self, diameter, wavelength, index, npts):

        self.diameter = diameter

        self.wavelength = wavelength

        self.index = index

        self.npts = npts

        self.Full = Meshes(ThetaBound=[0, 360],
                           PhiBound=[0, 360],
                           npts=npts)

        self.computeS1S2()

    def PlotFields(self):

        ax0, ax1, ax2 = self.GenFigure()

        S1 = self.GenS1Plot(ax0)

        S2 = self.GenS2Plot(ax1)

        SPF = self.Make3DField(ax2)

        plt.show()

        Plot = S1S2Plot(np.abs(self.S1), np.abs(self.S2), *SPF, self.Full)

    def GenS1Plot(self, axe):
        """Generate the polar subplot for S1.

        Parameters
        ----------
        axe : type
            subplot axe for S1.

        """

        data = (np.abs(self.S1))


        return data


    def GenS2Plot(self, axe):
        """Generate the polar subplot for S2.

        Parameters
        ----------
        axe : type
            subplot axe for S2.

        """
        data = (np.abs(self.S2))
        axe.plot(self.Full.PhiVec.Radian,
                 data,
                 'k')

        axe.fill_between(self.Full.PhiVec.Radian,
                         0,
                         data,
                         color='C1',
                         alpha=0.4)

        return data


    def Make3DField(self, axe):

        Theta = np.linspace(0, 2*np.pi, self.npts)

        Phi = np.linspace(0, np.pi, self.npts)

        FuncSPF = interp.interp2d(self.Field.Meshes.ThetaVec.Radian,
                                  self.Field.Meshes.PhiVec.Radian,
                                  self.Field.SPF,
                                  kind='cubic')

        SPF = FuncSPF(Theta, Phi)

        PHI, THETA = np.meshgrid(Theta, Phi)

        PHI, THETA = THETA, PHI

        X = SPF * np.sin(PHI) * np.cos(THETA)
        Y = SPF * np.sin(PHI) * np.sin(THETA)
        Z = SPF * np.cos(PHI)

        axe.plot_surface(Y,
                         X,
                         Z,
                         rstride=3,
                         cstride=3,
                         linewidth=0.5,
                         cmap=cm.bone,
                         antialiased=False,
                         alpha=1)

        axe.set_xlabel('X direction')
        axe.set_ylabel('Y direction')
        axe.set_zlabel('Z direction')
        norm = np.sqrt(np.max(X**2+Y**2))

        axe.set_xlim([-norm, norm])
        axe.set_ylim([-norm, norm])

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

            Parallel = np.outer(self.S1, np.sin(self.Full.PhiVec.Radian+Polarization))

            Perpendicular = np.outer(self.S2, np.cos(self.Full.PhiVec.Radian+Polarization))

        elif PolarizationAngle is None:

            self.Polarization = "None"

            Parallel = np.outer(self.S1,  np.ones(self.npts)/np.sqrt(2))

            Perpendicular = np.outer(self.S2, np.ones(self.npts)/np.sqrt(2))

        self.Field = Field(Perpendicular=Perpendicular,
                           Parallel=Parallel,
                           Meshes=self.Full
                           )

    def GenFigure(self):

        fig = plt.figure(figsize=(15, 5))

        ax0 = fig.add_subplot(131, projection='polar')
        ax1 = fig.add_subplot(132, projection='polar')
        ax2 = fig.add_subplot(133, projection='3d')

        ax0.set_title(r'$S_1$ [logscale]')
        ax1.set_title(r'$S_2$ [logscale]')
        ax2.set_title(r'3D Phase function Polarization: {0}'.format(self.Polarization))

        return ax0, ax1, ax2

    def SampleField(self, ThetaBound, PhiBound, npts=100, ThetaPol=0):

        FuncParallelReal = interp.interp2d(self.Field.Meshes.ThetaVec.Degree,
                                           self.Field.Meshes.PhiVec.Degree,
                                           np.real(self.Field.Parallel),
                                           kind='cubic')

        FuncParallelImage = interp.interp2d(self.Field.Meshes.ThetaVec.Degree,
                                            self.Field.Meshes.PhiVec.Degree,
                                            np.imag(self.Field.Parallel),
                                            kind='cubic')

        FuncPerpendicularReal = interp.interp2d(self.Field.Meshes.ThetaVec.Degree,
                                                self.Field.Meshes.PhiVec.Degree,
                                                np.real(self.Field.Perpendicular),
                                                kind='cubic')

        FuncPerpendicularImage = interp.interp2d(self.Field.Meshes.ThetaVec.Degree,
                                                 self.Field.Meshes.PhiVec.Degree,
                                                 np.imag(self.Field.Perpendicular),
                                                 kind='cubic')

        self.SampleMesh = Meshes(ThetaBound=ThetaBound,
                                 PhiBound=PhiBound,
                                 npts=npts)

        Perpendicular = FuncPerpendicularReal(self.SampleMesh.ThetaVec.DegreeMod, self.SampleMesh.PhiVec.DegreeMod)\
            + i * FuncPerpendicularImage(self.SampleMesh.ThetaVec.DegreeMod,
                                         self.SampleMesh.PhiVec.DegreeMod)

        Parallel = FuncParallelReal(self.SampleMesh.ThetaVec.DegreeMod, self.SampleMesh.PhiVec.DegreeMod)\
            + i * FuncParallelImage(self.SampleMesh.ThetaVec.DegreeMod,
                                    self.SampleMesh.PhiVec.DegreeMod)

        self.Sample = Field(Perpendicular=Perpendicular,
                            Parallel=Parallel,
                            Meshes=self.SampleMesh
                            )

    # -
