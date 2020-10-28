
from numpy import abs, arctan, angle, sqrt, real, imag, conjugate, exp, array, max
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patches as mpatches


class Field(object):

    def __init__(self,
                 Perpendicular: np.array,
                 Parallel: np.array,
                 Meshes: np.array):
        """
        Source -- https://www.physlab.org/wp-content/uploads/2016/07/Ch6-BYUOpticsBook_2013.pdf

        """
        self.__dict__ = Meshes.__dict__.copy()

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
            self.__SPF = self.ComputeSPF()
            return self.__SPF

        else:
            return self.__SPF


    @property
    def Stokes(self) -> None:
        if self.__Stokes is None:
            self.__Stokes = self.ComputeStokes()
            return self.__Stokes

        else:
            return self.__Stokes


    @property
    def Jones(self) -> None:
        if self.__Jones is None:
            self.__Jones = self.ComputeJones()
            return self.__Jones

        else:
            return self.__Jones


    def ComputeTotal(self) -> None:
        return sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)# * exp(complex(0,1)*self.delta)


    def ComputeDelay(self) -> None:
        return arctan(abs(self.Parallel)/abs(self.Perpendicular))


    def ComputeSPF(self) -> None:
        return abs(self.Parallel)**2 + abs(self.Perpendicular)**2


    def ComputeStokes(self) -> None:
        I = abs(self.Parallel)**2 + abs(self.Perpendicular)**2

        return array( [ I,
                        (abs(self.Parallel)**2 - abs(self.Perpendicular)**2)/I,
                        2*real( self.Parallel * conjugate(self.Perpendicular))/I,
                        -2*imag( self.Parallel * conjugate(self.Perpendicular))/I
                        ] )


    def ComputeJones(self) -> None:
        delta = angle(self.Parallel)-angle(self.Perpendicular)

        A = abs(self.Parallel) / sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        B = abs(self.Perpendicular) / sqrt(abs(self.Parallel)**2 + abs(self.Perpendicular)**2)

        return array([A, B*exp(complex(0,1)*delta)])


    def ComputeS1S2(self, CacheTrunk, Diameter, Wavelength) -> None:

        MuList = np.cos(self.Meshes.Phi.Vector.Radian)

        if CacheTrunk: MuList = np.round(MuList, CacheTrunk)

        self.S1, self.S2 = [], []

        for Mu in MuList:

            SizeParam = np.pi * Diameter / Wavelength

            S1, S2 = self.WrapS1S2(Mu)

            self.S1.append(S1)
            self.S2.append(S2)

        return S1, S2



    def GenStokesFig(self):
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        axes[0].set_title(r'$S_0$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$] ')

        axes[1].set_title(r'$S_3$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$]')

        [ ax.set_ylabel(r'Angle $\theta$') for ax in axes ]

        [ ax.set_xlabel(r'Angle $\phi$') for ax in axes ]

        return fig, axes


    def Plot(self, Function):
        pass


    def PlotStokes(self, RectangleTheta=[0,360], RectanglePhi=[0,360]):
        fig, axes = self.GenStokesFig()

        Origin = (RectangleTheta[1], RectanglePhi[1])

        SolidAngle = (RectangleTheta[0] - RectangleTheta[1], RectanglePhi[0] - RectanglePhi[1])

        n=8

        for i, ax in enumerate(axes):

            im = ax.pcolormesh(self.Phi.Vector.Degree,
                               self.Theta.Vector.Degree,
                               self.Stokes[3].T,
                               shading='auto')

            cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')

            cbar.ax.tick_params(labelsize='small')

            cbar.ax.locator_params(nbins=3)

            ax.quiver(self.Phi.Mesh.Degree[::n, ::n],
                      self.Theta.Mesh.Degree[::n, ::n],
                      self.Stokes[2][::n, ::n],
                      self.Stokes[1][::n, ::n],
                      units='width',
                      width=0.0005,
                      headwidth=30,
                      headlength=20,
                      headaxislength=20)

        [ax.add_patch(patches.Rectangle(Origin,*SolidAngle,linewidth=1,edgecolor='r',facecolor='none')) for ax in axes ]

        plt.show()














# -
