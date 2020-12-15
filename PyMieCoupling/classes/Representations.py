import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from typing import Tuple
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes
import matplotlib.ticker as tick
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from mpl_toolkits.axes_grid1 import make_axes_locatable
import polarTransform
from PyMieCoupling.functions.converts import NA2Angle


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

        ax0.contourf(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.real(self.Perpendicular))

        ax1.contourf(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.imag(self.Perpendicular))

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
    def __new__(cls, Field):

        cls.Meshes = Field.Meshes

        scamap = plt.cm.ScalarMappable(cmap='jet')

        cls.fcolors = scamap.to_rgba(np.real(Field.Perpendicular) + np.real(Field.Parallel) )

        SPF3D = cls.Make3D(cls,
                            Field.Parallel.__abs__()**2 + Field.Perpendicular.__abs__()**2,
                            cls.Meshes.Theta.Mesh.Radian,
                            cls.Meshes.Phi.Mesh.Radian,
                            )

        this = np.array(SPF3D, copy=False)

        this = np.asarray(this).view(cls)

        return this


    def __array_finalize__(self, obj):
        pass


    def __init__(self, Field):

        pass


    def Plot(cls):

        fig, ax = cls.GenFig()

        ax.plot_surface(*cls,
                         rstride     = 2,
                         cstride     = 2,
                         linewidth   = 0.003,
                         facecolors=cls.fcolors,
                         edgecolors='k',
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








class LP_FarField(object):

    def __init__(self,
                 Input,
                 Size,
                 Npts,
                 NA):

        self.Cartesian = Input
        self.Size = Size
        self.Npts = Npts
        self.X = np.linspace(-self.Size/2, self.Size/2, self.Npts)*1e6
        self.Y = np.linspace(-self.Size/2, self.Size/2, self.Npts)*1e6
        self._NA = NA
        self._PhiBound =  np.asarray( [0, NA2Angle(self.NA)] )
        self._ThetaBound = np.asarray([-180, 180])
        self.GetSpherical()
        self.Meshes = AngleMeshes(ThetaBound = self._ThetaBound,
                                  PhiBound   = self._PhiBound,
                                  ThetaNpts  = self.Spherical.shape[0],
                                  PhiNpts    = self.Spherical.shape[1])


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



    def Plot(self):
        fig = plt.figure(figsize=(12,3))
        ax0 = fig.add_subplot(121, projection = 'aitoff')
        ax1 = fig.add_subplot(122, projection = 'aitoff')

        ax0.contourf(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.real(self.Spherical.T))

        ax1.contourf(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.imag(self.Spherical.T))

        ax0.set_title('Real Part\n Far-Field Spherical Coordinates')
        ax0.set_xlabel(r'Angle $\phi$ [Degree]')
        ax0.set_ylabel(r'Angle $\theta$ [Degree]')
        ax0.grid()

        ax1.set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        ax1.set_xlabel(r'Angle $\phi$ [Degree]')
        ax1.set_ylabel(r'Angle $\theta$ [Degree]')
        ax1.grid()

        fig.tight_layout()



    def PlotPlan(self):
        fig = plt.figure(figsize=(12,3))
        ax0 = fig.add_subplot(141)
        ax1 = fig.add_subplot(142)
        ax2 = fig.add_subplot(143)
        ax3 = fig.add_subplot(144)

        ax0.pcolormesh(self.X,
                       self.Y,
                       self.Cartesian.real,
                       shading='auto')

        ax0.set_title('Real Part\n Far-Field Cartesian Coordinates')
        ax0.set_xlabel(r'X-Distance c * x [u.a.]')
        ax0.set_ylabel(r'Y-Distance c * y [u.a.]')

        ax1.pcolormesh(self.X,
                       self.Y,
                       self.Cartesian.imag,
                       shading='auto')

        ax1.set_title('Imaginary Part\n Far-Field Cartesian Coordinates')
        ax1.set_xlabel(r'X-Distance c * x [u.a.]')
        ax1.set_ylabel(r'Y-Distance c * y [u.a.]')

        ax2.pcolormesh(self.Meshes.Phi.Vector.Degree,
                       self.Meshes.Theta.Vector.Degree,
                       self.Spherical.real,
                       shading='auto')

        ax2.set_title('Real Part\n Far-Field Spherical Coordinates')
        ax2.set_xlabel(r'Angle $\phi$ [Degree]')
        ax2.set_ylabel(r'Angle $\theta$ [Degree]')

        ax3.pcolormesh(self.Meshes.Phi.Vector.Degree,
                       self.Meshes.Theta.Vector.Degree,
                       self.Spherical.imag,
                       shading='auto')

        ax3.set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        ax3.set_xlabel(r'Angle $\phi$ [Degree]')
        ax3.set_ylabel(r'Angle $\theta$ [Degree]')

        fig.tight_layout()


    @property
    def ThetaBound(self):
        return self._ThetaBound

    @property
    def PhiBound(self):
        return self._PhiBound

    @property
    def PhiOffset(self):
        return self.__PhiOffset

    @property
    def ThetaOffset(self):
        return self.__ThetaOffset

    @property
    def NA(self):
        return self._NA

    @NA.setter
    def NA(self, val: float):
        self._NA = val
        self.PhiBound =  np.asarray( [0, NA2Angle(self._NA)] )


    @ThetaBound.setter
    def ThetaBound(self, val: list):
        self._ThetaBound = np.asarray( val )
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = 0)
    @PhiBound.setter
    def PhiBound(self, val: list):
        self._PhiBound = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = 0)

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = val)

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self._ThetaOffset = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = val,
                                  PhiOffset          = 0)






class Detector_FarField(object):

    def __init__(self,
                 Npts,
                 NA,
                 ThetaOffset = 0,
                 PhiOffset   = 0):


        self.Npts = Npts
        self._NA = NA
        self._PhiBound =  np.asarray( [0, NA2Angle(self._NA)] )
        self._ThetaBound = np.asarray([-180, 180])
        self.__ThetaOffset, self.__PhiOffset = ThetaOffset, PhiOffset
        self.GetSpherical()
        self.Meshes = AngleMeshes(ThetaBound = self._ThetaBound + self.__ThetaOffset,
                                  PhiBound   = self._PhiBound + self.__PhiOffset,
                                  ThetaNpts  = self.Npts,
                                  PhiNpts    = self.Npts)


    def GetSpherical(self):
        self.Spherical = np.ones([self.Npts, self.Npts]) / (self.Npts*self.Npts)


    def Plot(self):
        fig = plt.figure(figsize=(8,4))
        ax0 = fig.add_subplot(121)
        ax1 = fig.add_subplot(122)

        im0 = ax0.pcolormesh(self.Meshes.Phi.Vector.Degree,
                             self.Meshes.Theta.Vector.Degree,
                             self.Spherical.real,
                             shading='auto')

        ax0.set_title('Real Part\n Far-Field Spherical Coordinates')
        ax0.set_xlabel(r'Angle $\phi$ [Degree]')
        ax0.set_ylabel(r'Angle $\theta$ [Degree]')

        im1 = ax1.pcolormesh(self.Meshes.Phi.Vector.Degree,
                             self.Meshes.Theta.Vector.Degree,
                             self.Spherical.imag,
                             shading='auto')

        ax1.set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        ax1.set_xlabel(r'Angle $\phi$ [Degree]')
        ax1.set_ylabel(r'Angle $\theta$ [Degree]')

        divider = make_axes_locatable(ax0)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im0, cax=cax, orientation='vertical')

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im1, cax=cax, orientation='vertical')

        fig.tight_layout()


    @property
    def NA(self):
        return self._NA

    @NA.setter
    def NA(self, val: float):
        self._NA = val
        self.PhiBound =  np.asarray( [0, NA2Angle(self._NA)] )


    @property
    def ThetaBound(self):
        return self._ThetaBound

    @property
    def PhiBound(self):
        return self._PhiBound

    @property
    def PhiOffset(self):
        return self.__PhiOffset

    @property
    def ThetaOffset(self):
        return self.__ThetaOffset

    @ThetaBound.setter
    def ThetaBound(self, val: list):
        self._ThetaBound = np.asarray( val )
        self.Meshes = AngleMeshes(ThetaNpts  = self.Npts,
                                  PhiNpts    = self.Npts,
                                  ThetaBound = self._ThetaBound,
                                  PhiBound   = self._PhiBound)

    @PhiBound.setter
    def PhiBound(self, val: list):
        self._PhiBound = val
        print(self._PhiBound)

        self.Meshes = AngleMeshes(ThetaNpts  = self.Npts,
                                  PhiNpts    = self.Npts,
                                  ThetaBound = self._ThetaBound,
                                  PhiBound   = self._PhiBound)



    @PhiOffset.setter
    def PhiOffset(self, val):
        self.__PhiOffset = val
        self.Meshes = AngleMeshes(ThetaNpts  = self.Npts,
                                  PhiNpts    = self.Npts,
                                  ThetaBound = self._ThetaBound,
                                  PhiBound   = self._PhiBound + val)


    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.__ThetaOffset = val
        self.Meshes = AngleMeshes(ThetaNpts  = self.Npts,
                                  PhiNpts    = self.Npts,
                                  ThetaBound = self._ThetaBound + val,
                                  PhiBound   = self._PhiBound)








# -
