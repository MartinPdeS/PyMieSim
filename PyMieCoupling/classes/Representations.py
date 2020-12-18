import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple
from PyMieCoupling.classes.Meshes import AngleMeshes

import matplotlib.ticker as tick
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

from PyMieCoupling.cpp.S1S2 import GetFields as Fields_CPP
from ai import cs


try:
    from PyMieCoupling.cpp.S1S2 import GetS1S2
except:
    try:
        from PyMieCoupling.cython.S1S2 import GetS1S2
    except:
        try:
            from PyMieCoupling.cython.S1S2 import GetS1S2
        except: ImportError


class Stokes(np.ndarray):
    def __new__(cls, Field):

        cls.Meshes = Field.Meshes

        Stokes = cls.GetStokes(cls, Field.Parallel, Field.Perpendicular)

        this = np.array(Stokes, copy=False)

        this = np.asarray(this).view(cls)

        return this


    def __array_finalize__(self, obj):
        pass


    def __init__(self, Field):
        pass


    def GenFig(self):

        fig, axes = plt.subplots(1, 2, figsize=(6, 3))

        axes[0].set_title('$S_0$ Stokes parameter of Far-Field \n Projection on [$S_1, S_2$] ')

        axes[1].set_title('$S_3$ Stokes parameter of Far-Field \n Projection on [$S_1, S_2$]')

        [ ax.set_ylabel(r'Angle $\theta$ [Degree]') for ax in axes ]

        [ ax.set_xlabel(r'Angle $\phi$ [Degree]') for ax in axes ]

        return fig, axes


    def Plot(cls):

        fig, axes = cls.GenFig()
        n = 6
        for n, ax in enumerate(axes):
            m = n * 3

            im = ax.pcolormesh(cls.Meshes.Phi.Vector.Degree,
                               cls.Meshes.Theta.Vector.Degree,
                               cls.Array[m,:,:],
                               shading='auto',
                               )

            ax.streamplot(cls.Meshes.Phi.Mesh.Degree[::n, ::n].T,
                          cls.Meshes.Theta.Mesh.Degree[::n, ::n].T,
                          cls.Array[1,::n,::n],
                          cls.Array[2,::n,::n])

            cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')
            cbar.ax.tick_params(labelsize='small')
            cbar.ax.locator_params(nbins=3)

        plt.show(block=False)


    def GetStokes(cls, Parallel, Perpendicular):

        Array = np.empty( [4, *Parallel.shape] )

        I = Parallel.__abs__()**2 + Perpendicular.__abs__()**2
        Array[0,:,:] = I
        Array[1,:,:] = (Parallel.__abs__()**2 - Perpendicular.__abs__()**2)/I
        Array[2,:,:] = 2 * ( Parallel * Perpendicular.conjugate() ).real / I
        Array[3,:,:] = -2 * ( Parallel.conjugate() * Perpendicular ).imag / I

        cls.Array = Array



class Field(object):

    def __init__(self,
                 Perpendicular: np.ndarray,
                 Parallel:      np.ndarray,
                 Meshes:        AngleMeshes):
        """
        Source -- https://www.physlab.org/wp-content/uploads/2016/07/Ch6-BYUOpticsBook_2013.pdf

        """
        self.__dict__ = Meshes.__dict__
        self.Perpendicular, self.Parallel = Perpendicular, Parallel
        self.Meshes = Meshes


    def Plot(self):
        fig, axes = plt.subplots(nrows      = 1,
                                 ncols      = 2,
                                 figsize    = (12,3),
                                 subplot_kw = {'projection':'mollweide'})

        axes[0].pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian-np.pi/2,
                     np.real(self.Perpendicular),
                     shading='auto')

        axes[1].pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian-np.pi/2,
                     np.imag(self.Perpendicular),
                     shading='auto')

        [ax.set_ylabel(r'Angle $\phi$ [Degree]') for ax in axes]
        [ax.set_xlabel(r'Angle $\theta$ [Degree]') for ax in axes]
        [ax.grid() for ax in axes]
        axes[0].set_title('Real Part\n Far-Field Spherical Coordinates')
        axes[1].set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        fig.tight_layout()




class SPF(np.ndarray):
    def __new__(cls, Index, SizeParam, Polarization=0, num=201):
        thetaVec = np.linspace(-np.pi, np.pi, num)
        phiVec   = np.linspace(0, np.pi, num)

        thetaMesh, phiMesh = np.meshgrid(thetaVec, phiVec)

        Parallel, Perpendicular = Fields_CPP(Index,
                                             SizeParam,
                                             thetaMesh.flatten(),
                                             phiMesh.flatten(),
                                             phiVec,
                                             num,
                                             num,
                                             Polarization  = Polarization);

        scamap = plt.cm.ScalarMappable(cmap='jet')
        cls.fcolors = scamap.to_rgba(np.real(Perpendicular) + np.real(Parallel) )

        SPF3D = cs.sp2cart(Parallel.__abs__()**2 + Perpendicular.__abs__()**2,
                           phiMesh - np.pi/2,
                           thetaMesh,
                           )

        this = np.array(SPF3D, copy=False)
        this = np.asarray(this).view(cls)

        return this


    def Plot(cls):

        fig, ax = cls.GenFig()

        ax.plot_surface(*cls,
                         rstride     = 2,
                         cstride     = 2,
                         linewidth   = 0.003,
                         facecolors  = cls.fcolors,
                         edgecolors  = 'k',
                         antialiased = False,
                         )

        xLim = ax.get_xlim(); yLim = ax.get_ylim(); zLim = ax.get_zlim()

        Min = min(xLim[0],yLim[0]); Max = max(xLim[1], yLim[1])

        ax.set_xlim(Min, Max)

        ax.set_ylim(Min, Max)

        plt.show(block=False)


    def GenFig(self) -> Tuple[plt.figure, plt.axes]:

        fig, ax = plt.subplots(1, figsize=(3, 3), subplot_kw = {'projection':'3d'})
        ax.set_title(r'Complex Scattering Phase Function: Real{$ E_{||}$}')
        ax.set_ylabel(r'Y-direction')
        ax.set_xlabel(r'X-direction')
        ax.set_zlabel(r'Z-direction')

        return fig, ax



class S1S2(np.ndarray):
    """Short summary.

    Parameters
    ----------
    SizeParam : np.array
        Description of parameter `SizeParam`.
    Index : float
        Description of parameter `Index`.
    Meshes : AngleMeshes
        Description of parameter `Meshes`.
    CacheTrunk : bool
        Description of parameter `CacheTrunk`.

    Attributes
    ----------
    Array : type
        Description of attribute `Array`.
    Meshes
    SizeParam
    Index

    """

    def __new__(cls,
                SizeParam:  np.array,
                Index:      float,
                Meshes:     AngleMeshes):

        temp = GetS1S2(Index,
                       SizeParam,
                       Meshes.Phi.Vector.Radian)

        this = np.array(temp, copy=False)
        this = np.asarray(this).view(cls)
        return this


    def __init__(self, SizeParam, Index, Meshes):
        self.Meshes = Meshes


    def Plot(self) -> None:

        fig, axes = plt.subplots(nrows      = 1,
                                 ncols      = 2,
                                 subplot_kw = {'projection':'polar'})

        axes[0].set_title('S1 function'); axes[1].set_title('S2 function')

        data = np.abs(self)

        for ni, ax in enumerate(axes):

            ax.plot(self.Meshes.Phi.Vector.Radian,
                    data[ni],
                    'k')

            ax.fill_between(self.Meshes.Phi.Vector.Radian,
                            0,
                            data[ni],
                            color='C0',
                            alpha=0.4)












# -
