import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as patches
import matplotlib.patches as mpatches
from PyMieCoupling.functions.Misc import Make3D
import functools
from typing import Tuple
import PyMieScatt


class Stokes(object):
    def __init__(self,
                 Parallel: np.array,
                 Perpendicular: np.array,
                 Meshes) -> None:

        self.shape = np.shape(Parallel)

        self.Meshes = Meshes

        self.Array = np.empty( [4, *self.shape] )

        I = abs(Parallel)**2 + abs(Perpendicular)**2

        self.Array[0,:,:] = I

        self.Array[1,:,:] = (abs(Parallel)**2 - abs(Perpendicular)**2)/I

        self.Array[2,:,:] = 2*np.real( Parallel * np.conjugate(Perpendicular))/I

        self.Array[3,:,:] = -2*np.imag( Parallel * np.conjugate(Perpendicular))/I


    def GenFig(self):

        fig, axes = plt.subplots(1, 2, figsize=(6, 3))

        axes[0].set_title('$S_0$ Stokes parameter of Far-Field \n Projection on [$S_1, S_2$] ')

        axes[1].set_title('$S_3$ Stokes parameter of Far-Field \n Projection on [$S_1, S_2$]')

        [ ax.set_ylabel(r'Angle $\theta$') for ax in axes ]

        [ ax.set_xlabel(r'Angle $\phi$') for ax in axes ]

        return fig, axes



    def Plot(self,
             RectangleTheta = [-5,5],
             RectanglePhi = [-5,5]):

        fig, axes = self.GenFig()

        Origin = (RectangleTheta[1], RectanglePhi[1])

        SolidAngle = (RectangleTheta[0] - RectangleTheta[1], RectanglePhi[0] - RectanglePhi[1])

        n=8

        for i, ax in enumerate(axes):

            im = ax.pcolormesh(self.Meshes.Phi.Vector.Degree,
                               self.Meshes.Theta.Vector.Degree,
                               self.Array[3].T,
                               shading='auto')

            cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')

            cbar.ax.tick_params(labelsize='small')

            cbar.ax.locator_params(nbins=3)

            ax.quiver(self.Meshes.Phi.Mesh.Degree[::n, ::n],
                      self.Meshes.Theta.Mesh.Degree[::n, ::n],
                      self.Array[2][::n, ::n],
                      self.Array[1][::n, ::n],
                      units='width',
                      width=0.0005,
                      headwidth=30,
                      headlength=20,
                      headaxislength=20)

        [ax.add_patch(patches.Rectangle(Origin,*SolidAngle,linewidth=1,edgecolor='r',facecolor='none')) for ax in axes ]

        plt.show()



    def __repr__(self):

        s0, s1 = self.shape

        return '\nStokes Field representation    \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(s0,
                        s1,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )



class Jones(object):
    def __init__(self,
                 Parallel: np.array,
                 Perpendicular: np.array,
                 Meshes) -> None:

        self.shape = np.shape(Parallel)

        self.Meshes = Meshes

        self.Array = np.empty( [2, *self.shape] )

        delta = np.angle(Parallel)-np.angle(Perpendicular)

        A = np.abs(Parallel) / np.sqrt(abs(Parallel)**2 + np.abs(Perpendicular)**2)

        B = np.abs(Perpendicular) / np.sqrt(abs(Parallel)**2 + np.abs(Perpendicular)**2)

        self.Array = np.array([A, B*np.exp(complex(0,1)*delta)])



    def GenFig(self):

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        axes[0].set_title(r'$S_0$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$] ')

        axes[1].set_title(r'$S_3$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$]')

        [ ax.set_ylabel(r'Angle $\theta$') for ax in axes ]

        [ ax.set_xlabel(r'Angle $\phi$') for ax in axes ]

        return fig, axes



    def Plot(self,
             RectangleTheta = [-5,5],
             RectanglePhi = [-5,5]):

        fig, axes = self.GenFig()

        Origin = (RectangleTheta[1], RectanglePhi[1])

        SolidAngle = (RectangleTheta[0] - RectangleTheta[1], RectanglePhi[0] - RectanglePhi[1])





    def __repr__(self):

        s0, s1 = self.shape

        return '\Jones Field representation    \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(s0,
                        s1,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )








class SPF(object):
    def __init__(self,
                 Parallel: np.array,
                 Perpendicular: np.array,
                 Meshes) -> None:

        self.shape = np.shape(Parallel)

        self.Meshes = Meshes

        self.Array = np.abs(Parallel)**2 + np.abs(Perpendicular)**2



    def GenFig(self):

        fig = plt.figure(figsize=(3, 3))

        ax = fig.add_subplot(111, projection='3d')

        ax.set_title(r'Scattering Phase Function')

        ax.set_ylabel(r'Y-direction')

        ax.set_xlabel(r'X-direction')

        ax.set_zlabel(r'Z-direction')

        return fig, ax



    def Plot(self,
             RectangleTheta = [-5,5],
             RectanglePhi = [-5,5]):

        fig, ax = self.GenFig()

        Origin = (RectangleTheta[1], RectanglePhi[1])

        SolidAngle = (RectangleTheta[0] - RectangleTheta[1], RectanglePhi[0] - RectanglePhi[1])

        SPF3D = Make3D(self.Array,
                       self.Meshes.Phi.Mesh.Radian,
                       self.Meshes.Theta.Mesh.Radian)

        ax.plot_surface(*SPF3D,
                         rstride=3,
                         cstride=3,
                         linewidth=0.5,
                         cmap=cm.bone,
                         antialiased=False,
                         alpha=1)

        plt.show()


    def __repr__(self):

        s0, s1 = self.shape

        return '\nScattering Phase Function      \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(s0,
                        s1,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )







class S1S2(object):
    def __init__(self,
                 SizeParam: np.array,
                 MuList: np.array,
                 Index: float,
                 Meshes,
                 CacheTrunk: bool = None) -> None:

        self.shape = np.shape(MuList)

        self.Meshes, self.SizeParam = Meshes, SizeParam

        self.Index = Index

        if CacheTrunk: MuList = np.round(MuList, CacheTrunk)

        S1, S2 = [], []

        for Mu in MuList:

            temp0, temp1 = self.WrapS1S2(Mu)

            S1.append(temp0)
            S2.append(temp1)

        self.Array = [S1, S2]


    @functools.lru_cache(maxsize=201)
    def WrapS1S2(self, Mu) -> Tuple[float, float]:

        S1, S2 = PyMieScatt.MieS1S2(self.Index,
                                    self.SizeParam,
                                    Mu)

        return S1, S2


    def GenFig(self):

        fig = plt.figure(figsize=(6, 3))

        ax0 = fig.add_subplot(121, projection = 'polar')

        ax1 = fig.add_subplot(122, projection = 'polar')

        ax0.set_title(r'S1 function')

        ax1.set_title(r'S2 function')

        return fig, [ax0, ax1]


    def Plot(self):

        fig, axes = self.GenFig()

        axes[0].plot(self.Meshes.Phi.Vector.Radian,
                     np.abs(self.Array[0]),
                     'k')

        axes[0].fill_between(self.Meshes.Phi.Vector.Radian,
                             0,
                             np.abs(self.Array[0]),
                             color='C0',
                             alpha=0.4)

        axes[1].plot(self.Meshes.Phi.Vector.Radian,
                     np.abs(self.Array[1]),
                     'k')

        axes[1].fill_between(self.Meshes.Phi.Vector.Radian,
                             0,
                             np.abs(self.Array[1]),
                             color='C1',
                             alpha=0.4)

        plt.show()


    def __repr__(self):

        s0, s1 = self.shape

        return '\nScattering Phase Function      \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(s0,
                        s1,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )





# -
