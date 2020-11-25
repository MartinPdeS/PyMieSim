
import numpy as np
from PyMieCoupling.classes.Meshes import ScatMeshes
from PyMieCoupling.classes.Misc import Source
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple




#from PyMieCoupling.cython.S1S2 import MieS1S2
#from PyMieCoupling.python.S1S2 import MieS1S2# as S1S2_PYTHON
from PyMieCoupling.cpp.S1S2 import MieS1S2# as S1S2_CPP


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
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Npts:        int         = None,
                 Meshes:      ScatMeshes  = None,
                 ThetaBound:  list        = [-180, 180],
                 ThetaOffset: float       = 0,
                 PhiBound:    list        = [-180, 180],
                 PhiOffset:   float       = 0) -> None:


        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter / 2 )

        if Meshes:
            self.Meshes = Meshes
        else:
            self.GenMesh(ThetaBound, PhiBound, ThetaOffset, PhiOffset, Npts)

        self.__S1S2, self.__Field = None, None



    def GenMesh(self,
                ThetaBound:  list       = [-90, 90],
                PhiBound:    list       = [-90, 90],
                ThetaOffset: float      = 0,
                PhiOffset:   float      = 0,
                Npts:        int        = 101):



        self.Meshes = ScatMeshes(ThetaBound = np.array(ThetaBound) + ThetaOffset,
                                 PhiBound   = np.array(PhiBound) + PhiOffset,
                                 Npts       = Npts)


    @property
    def S1S2(self) -> np.ndarray:
        if self.__S1S2 is None:
            self.__S1S2 = RepS1S2(SizeParam  = self.SizeParam,
                                  Index      = self.Index,
                                  Meshes     = self.Meshes)
            return self.__S1S2

        else:
            return self.__S1S2


    @property
    def Field(self) -> ScatMeshes:
        if self.__Field is None:
            self.GenField()
            return self.__Field
        else:
            return self.__Field



    def GenField(self):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.

        """

        S1S2 = MieS1S2(self.Index,
                         self.SizeParam,
                         self.Meshes.Phi.Vector.Radian,
                         )

        Parallel = np.outer(S1S2[0], np.sin(self.Meshes.Theta.Vector.Radian))

        Perpendicular = np.outer(S1S2[1], np.cos(self.Meshes.Theta.Vector.Radian))

        self.__Field = Field(Perpendicular = Perpendicular,
                             Parallel      = Parallel,
                             Meshes        = self.Meshes
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






class Stokes(object):
    """Short summary.

    Parameters
    ----------
    Parallel : np.ndarray
        Description of parameter `Parallel`.
    Perpendicular : np.ndarray
        Description of parameter `Perpendicular`.
    Meshes : ScatMeshes
        Meshes of the scatterer.

    Attributes
    ----------
    Array : type
        Description of attribute `Array`.
    Meshes

    """
    def __init__(self,
                 Parallel:      np.ndarray,
                 Perpendicular: np.ndarray,
                 Meshes:        ScatMeshes) -> None:

        self.Meshes = Meshes

        self.Array = GetStokes(Parallel      = Parallel,
                               Perpendicular = Perpendicular)


    def __repr__(self) -> str:

        return '\nStokes Field representation    \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(*self.Meshes.Theta.Boundary.Degree.shape,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )


    def GenFig(self):

        fig, axes = plt.subplots(1, 2, figsize=(6, 3))

        axes[0].set_title('$S_0$ Stokes parameter of Far-Field \n Projection on [$S_1, S_2$] ')

        axes[1].set_title('$S_3$ Stokes parameter of Far-Field \n Projection on [$S_1, S_2$]')

        [ ax.set_ylabel(r'Angle $\theta$') for ax in axes ]

        [ ax.set_xlabel(r'Angle $\phi$') for ax in axes ]

        return fig, axes


    def Plot(self):

        fig, axes = self.GenFig()


        n = 3

        ax = axes[0]
        im = ax.pcolormesh(self.Meshes.Phi.Mesh.Degree,
                           self.Meshes.Theta.Mesh.Degree,
                           self.Array[0,:,:],
                           shading='auto',
                           )

        cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        ax.quiver(self.Meshes.Phi.Mesh.Degree[::n, ::n],
                  self.Meshes.Theta.Mesh.Degree[::n, ::n],
                  self.Array[1,::n,::n],
                  self.Array[2,::n,::n],
                  units          = 'width',
                  width          = 0.0005,
                  headwidth      = 30,
                  headlength     = 20,
                  headaxislength = 20)

        ax = axes[1]
        im = ax.pcolormesh(self.Meshes.Phi.Mesh.Degree,
                           self.Meshes.Theta.Mesh.Degree,
                           self.Array[3,:,:],
                           shading='auto',
                           )

        cbar = plt.colorbar(im, ax=ax, pad=0.15, orientation='horizontal')

        cbar.ax.tick_params(labelsize='small')

        cbar.ax.locator_params(nbins=3)

        ax.quiver( self.Meshes.Phi.Mesh.Degree[::n, ::n],
                  self.Meshes.Theta.Mesh.Degree[::n, ::n],
                  self.Array[1,::n,::n],
                  self.Array[2,::n,::n],
                  units          = 'width',
                  width          = 0.0005,
                  headwidth      = 30,
                  headlength     = 20,
                  headaxislength = 20)



        plt.show()


class Jones(object):
    """Short summary.

    Parameters
    ----------
    Parallel : np.ndarray
        Description of parameter `Parallel`.
    Perpendicular : np.ndarray
        Description of parameter `Perpendicular`.
    Meshes : ScatMeshes
        Meshes of the scatterer.

    Attributes
    ----------
    Array : type
        Description of attribute `Array`.
    Meshes

    """

    def __init__(self,
                 Parallel:      np.ndarray,
                 Perpendicular: np.ndarray,
                 Meshes:        ScatMeshes) -> None:

        self.Meshes = Meshes

        self.Array = ComputeJones(Parallel, Perpendicular)


    def __repr__(self) -> str:

        return '\nJones Field representation     \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(*self.Meshes.Theta.Boundary.Degree.shape,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )


    def GenFig(self) -> Tuple[plt.figure, plt.axes]:

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        axes[0].set_title(r'$S_0$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$] ')

        axes[1].set_title(r'$S_3$ Stokes parameter of Far-Field & Projection on [$S_1, S_2$]')

        [ ax.set_ylabel(r'Angle $\theta$') for ax in axes ]

        [ ax.set_xlabel(r'Angle $\phi$') for ax in axes ]

        return fig, axes




class SPF(object):
    """Short summary.

    Parameters
    ----------
    Parallel : np.ndarray
        Description of parameter `Parallel`.
    Perpendicular : np.ndarray
        Description of parameter `Perpendicular`.
    Meshes : ScatMeshes
        Meshes of the scatterer.

    Attributes
    ----------
    Array : type
        Description of attribute `Array`.
    Meshes

    """

    def __init__(self,
                 Parallel:      np.ndarray,
                 Perpendicular: np.ndarray,
                 Meshes:        ScatMeshes) -> None:

        self.Meshes = Meshes

        self.Array = Parallel.__abs__()**2 + Perpendicular.__abs__()**2


    def __repr__(self) -> str:

        return '\nScattering Phase Function      \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(*self.Meshes.Theta.Boundary.Degree.shape,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )


    def GenFig(self) -> Tuple[plt.figure, plt.axes]:

        fig = plt.figure(figsize=(3, 3))

        ax = fig.add_subplot(111, projection='3d')

        ax.set_title(r'Scattering Phase Function')

        ax.set_ylabel(r'Y-direction')

        ax.set_xlabel(r'X-direction')

        ax.set_zlabel(r'Z-direction')

        return fig, ax



    def Plot(self):

        fig, ax = self.GenFig()

        SPF3D = Make3D(self.Array,
                       self.Meshes.Phi.Mesh.Radian,
                       self.Meshes.Theta.Mesh.Radian)

        ax.plot_surface(*SPF3D,
                         rstride     = 3,
                         cstride     = 3,
                         linewidth   = 0.5,
                         cmap        = cm.bone,
                         antialiased = False,
                         alpha       = 1)

        plt.show()



class RepS1S2(object):
    """Short summary.

    Parameters
    ----------
    SizeParam : np.array
        Description of parameter `SizeParam`.
    Index : float
        Description of parameter `Index`.
    Meshes : ScatMeshes
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

    def __init__(self,
                 SizeParam:  np.array,
                 Index:      float,
                 Meshes:     ScatMeshes,
                 CacheTrunk: bool = None) -> None:

        self.Meshes, self.SizeParam = Meshes, SizeParam

        self.Index = Index

        self.S1S2 = MieS1S2(self.Index,
                            self.SizeParam,
                            self.Meshes.Phi.Vector.Radian,
                            )


    def GenFig(self) -> Tuple[plt.figure, plt.axes]:

        fig = plt.figure(figsize=(6, 3))

        ax0 = fig.add_subplot(121, projection = 'polar')

        ax1 = fig.add_subplot(122, projection = 'polar')

        ax0.set_title(r'S1 function')

        ax1.set_title(r'S2 function')

        return fig, [ax0, ax1]


    def Plot(self) -> None:

        fig, axes = self.GenFig()
        data = np.abs( self.S1S2 )

        for ni, ax in enumerate(axes):

            ax.plot(self.Meshes.Phi.Vector.Radian,
                    data[ni],
                    'k')

            ax.fill_between(self.Meshes.Phi.Vector.Radian,
                            0,
                            data[ni],
                            color='C0',
                            alpha=0.4)

        plt.show()


    def __repr__(self) -> str:

        return '\nScattering Phase Function      \
                \nField dimensions: {0}x{1}      \
                \nTheta boundary: {2} deg.       \
                \nPhi boundary: {3} deg'         \
                .format(*self.Meshes.Phi.Vector.Radian.shape,
                        self.Meshes.Theta.Boundary.Degree,
                        self.Meshes.Phi.Boundary.Degree )



def GetJones(Parallel:      np.ndarray,
             Perpendicular: np.ndarray) -> np.ndarray:

    Array = np.empty( [2, * Parallel.shape] )

    delta = np.angle(Parallel) - np.angle(Perpendicular)

    A = Parallel.__abs__() / np.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    B = Perpendicular.__abs__() / np.sqrt(Parallel.__abs__()**2 + Perpendicular.__abs__()**2)

    return np.array([A, B * np.exp(complex(0,1)*delta)])




def GetStokes(Parallel:      np.ndarray,
              Perpendicular: np.ndarray) -> np.ndarray:

    Array = np.empty( [4, *Parallel.shape] )

    I = Parallel.__abs__()**2 + Perpendicular.__abs__()**2
    Array[0,:,:] = I

    Array[1,:,:] = (Parallel.__abs__()**2 - Perpendicular.__abs__()**2)/I

    Array[2,:,:] = 2 * ( Parallel * Perpendicular.conjugate() ).real / I

    Array[3,:,:] = -2 * ( Parallel.conjugate() * Perpendicular ).imag / I

    return Array



def Make3D(item:      np.array,
           PhiMesh:   np.array,
           ThetaMesh: np.array) -> Tuple[np.array, np.array, np.array]:

    X = item * np.sin(PhiMesh) * np.cos(ThetaMesh)

    Y = item * np.sin(PhiMesh) * np.sin(ThetaMesh)

    Z = item * np.cos(PhiMesh)

    return X, Y, Z










class Field(object):

    def __init__(self,
                 Perpendicular: np.ndarray,
                 Parallel:      np.ndarray,
                 Meshes:        ScatMeshes):
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
            self.__SPF = SPF(Parallel      = self.Parallel,
                             Perpendicular = self.Perpendicular,
                             Meshes        = self.Meshes)
            return self.__SPF

        else:
            return self.__SPF


    @property
    def Stokes(self) -> None:
        if self.__Stokes is None:
            self.__Stokes = Stokes(Parallel      = self.Parallel,
                                   Perpendicular = self.Perpendicular,
                                   Meshes        = self.Meshes)
            return self.__Stokes

        else:
            return self.__Stokes


    @property
    def Jones(self) -> None:
        if self.__Jones is None:
            self.__Jones = Jones(Parallel      = self.Parallel,
                                 Perpendicular = self.Perpendicular,
                                 Meshes        = self.Meshes)
            return self.__Jones

        else:
            return self.__Jones


    def ComputeTotal(self) -> None:
        return np.sqrt(self.Parallel.__abs__()**2 +\
                    self.Perpendicular.__abs__()**2)   # * exp(complex(0,1)*self.delta)


    def ComputeDelay(self) -> None:
        return np.arctan( self.Parallel.__abs__() / self.Perpendicular.__abs__() )


    def ComputeSPF(self) -> None:
        return self.Parallel.__abs__()**2 + self.Perpendicular.__abs__()**2







# -




    # -
