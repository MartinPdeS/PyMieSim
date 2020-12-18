import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from typing import Tuple
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
from PyMieCoupling.classes.BaseClasses import BaseFarField
import matplotlib.ticker as tick
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from mpl_toolkits.axes_grid1 import make_axes_locatable
import polarTransform
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.cpp.S1S2 import GetFields as Fields_CPP

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

        ax = axes[0]
        im = ax.pcolormesh(cls.Meshes.Phi.Vector.Degree,
                           cls.Meshes.Theta.Vector.Degree,
                           cls.Array[0,:,:],
                           shading='auto',
                           )

        cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        ax.streamplot(cls.Meshes.Phi.Mesh.Degree[::n, ::n].T,
                      cls.Meshes.Theta.Mesh.Degree[::n, ::n].T,
                      cls.Array[1,::n,::n],
                      cls.Array[2,::n,::n],
                  )

        ax = axes[1]
        im = ax.pcolormesh(cls.Meshes.Phi.Vector.Degree,
                           cls.Meshes.Theta.Vector.Degree,
                           cls.Array[3,:,:],
                           shading='auto',
                           )

        cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        ax.streamplot(cls.Meshes.Phi.Mesh.Degree[::n, ::n].T,
                      cls.Meshes.Theta.Mesh.Degree[::n, ::n].T,
                      cls.Array[1,::n,::n],
                      cls.Array[2,::n,::n],
                  )

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

        fig = plt.figure(figsize=(12,3))
        ax0 = fig.add_subplot(121, projection = 'mollweide')
        ax1 = fig.add_subplot(122, projection = 'mollweide')

        ax0.pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian-np.pi/2,
                     np.real(self.Perpendicular),
                     shading='auto')

        ax1.pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian-np.pi/2,
                     np.imag(self.Perpendicular),
                     shading='auto')

        ax0.set_title('Real Part\n Far-Field Spherical Coordinates')
        ax0.set_ylabel(r'Angle $\phi$ [Degree]')
        ax0.set_xlabel(r'Angle $\theta$ [Degree]')
        ax0.grid()

        ax1.set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        ax1.set_ylabel(r'Angle $\phi$ [Degree]')
        ax1.set_xlabel(r'Angle $\theta$ [Degree]')
        ax1.grid()

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

        SPF3D = cls.Make3D(cls,
                            Parallel.__abs__()**2 + Perpendicular.__abs__()**2,
                            thetaMesh,
                            phiMesh,
                            )

        this = np.array(SPF3D, copy=False)

        this = np.asarray(this).view(cls)

        return this


    def __array_finalize__(self, obj):
        pass


    def __init__(self, Index, SizeParam, Polarization=0):

        pass


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

        ax.set_xlim(Min, Max )

        ax.set_ylim(Min, Max )

        plt.show(block=False)


    def Make3D(cls,
               item:      np.array,
               ThetaMesh: np.array,
               PhiMesh:   np.array,
               ) -> Tuple[np.array, np.array, np.array]:

        PhiMesh = PhiMesh %np.pi

        X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

        Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

        Z = item * np.cos(PhiMesh)

        return X, Y, Z


    def GenFig(self) -> Tuple[plt.figure, plt.axes]:

        fig = plt.figure(figsize=(3, 3))

        ax = fig.add_subplot(111, projection='3d')

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


    def __array_finalize__(self, obj):
        pass


    def __init__(self, SizeParam, Index, Meshes):

        self.Meshes = Meshes



    def GenFig(self) -> Tuple[plt.figure, plt.axes]:

        fig = plt.figure(figsize=(6, 3))

        ax0 = fig.add_subplot(121, projection = 'polar')

        ax1 = fig.add_subplot(122, projection = 'polar')

        ax0.set_title(r'S1 function')

        ax1.set_title(r'S2 function')

        return fig, [ax0, ax1]



    def Plot(self) -> None:

        fig, axes = self.GenFig()

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

        plt.show(block=False)









class LP_FarField(BaseFarField):

    def __init__(self,
                 Input,
                 Size:        float,
                 Npts:        int     = 101,
                 NA:          float   = 0.2,
                 PhiOffset:   float   = 0,
                 ThetaOffset: float   = 0):

        self.Cartesian = Input
        self.Size, self.Npts, self._NA  = Size, Npts, NA
        self._PhiBound, self._ThetaBound  =  np.asarray( [0, NA2Angle(self._NA)] ), np.asarray([-180, 180])
        self._PhiOffset, self._ThetaOffset = PhiOffset, ThetaOffset
        self.GetSpherical()

    def GetSpherical(self):
        polarImageReal, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.real,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = 40,
                                                                        finalAngle    = 2*np.pi)

        polarImageimag, ptSettings = polarTransform.convertToPolarImage(self.Cartesian.imag,
                                                                        center        = [self.Npts//2, self.Npts//2],
                                                                        initialRadius = 0,
                                                                        finalRadius   = 40,
                                                                        finalAngle    = 2*np.pi)

        self.Spherical = polarImageReal + complex(0,1) * polarImageimag

        self.Meshes = AngleMeshes(ThetaBound  = self._ThetaBound,
                                  PhiBound    = self._PhiBound,
                                  ThetaNpts   = self.Spherical.shape[0],
                                  PhiNpts     = self.Spherical.shape[1],
                                  PhiOffset   = self._PhiOffset,
                                  ThetaOffset = self._ThetaOffset)



class Detector_FarField(BaseFarField):

    def __init__(self,
                 Npts:        int   = 101,
                 NA:          float = 0.2,
                 ThetaOffset: float = 0,
                 PhiOffset:   float = 0):


        self.Npts, self._NA = Npts, NA

        self._PhiBound, self._ThetaBound  =  np.asarray( [0, NA2Angle(self._NA)] ), np.asarray([-180, 180])

        self._PhiOffset, self._ThetaOffset = PhiOffset, ThetaOffset

        self.GetSpherical()



    def GetSpherical(self):
        self.Spherical = np.ones([self.Npts, self.Npts]) / (self.Npts*self.Npts)

        self.Meshes = AngleMeshes(ThetaBound  = self._ThetaBound,
                                  PhiBound    = self._PhiBound,
                                  ThetaNpts   = self.Npts,
                                  PhiNpts     = self.Npts,
                                  PhiOffset   = self._PhiOffset,
                                  ThetaOffset = self._ThetaOffset)






# -
