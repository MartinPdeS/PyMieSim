import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.cpp.S1S2 import GetFieldsFromMesh
from PyMieCoupling.utils import PlotStructuredSphere


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

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

            im = ax.pcolormesh(cls.Meshes.Phi.Degree,
                               cls.Meshes.Theta.Degree,
                               cls.Array[m,:,:],
                               shading='auto',
                               )

            ax.streamplot(cls.Meshes.Phi.Degree[::n, ::n].T,
                          cls.Meshes.Theta.Degree[::n, ::n].T,
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
                 Parent):
        """
        Source -- https://www.physlab.org/wp-content/uploads/2016/07/Ch6-BYUOpticsBook_2013.pdf

        """
        self.Perpendicular, self.Parallel = Perpendicular, Parallel
        self.Parent = Parent


    def Plot(self):
        fig, axes = plt.subplots(nrows      = 1,
                                 ncols      = 2,
                                 figsize    = (12,3),
                                 subplot_kw = {'projection':'mollweide'})

        axes[0].pcolormesh(
                         self.Parent.Meshes.Theta.Radian,
                         -( self.Parent.Meshes.Phi.Radian - np.pi/2),
                         np.real(self.Perpendicular),
                         shading='auto')

        axes[1].pcolormesh(
                         self.Parent.Meshes.Theta.Radian,
                         -( self.Parent.Meshes.Phi.Radian - np.pi/2),
                         np.imag(self.Perpendicular),
                         shading='auto')

        [ax.set_ylabel(r'Angle $\phi$ [Degree]') for ax in axes]
        [ax.set_xlabel(r'Angle $\theta$ [Degree]') for ax in axes]
        [ax.grid() for ax in axes]
        axes[0].set_title('Real Part\n Far-Field Spherical Coordinates')
        axes[1].set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        fig.tight_layout()




class SPF(dict):

    def __init__(self, Parent, Num):

        self['Phi'], self['Theta'] = np.mgrid[-np.pi/2:np.pi/2:complex(Num), -np.pi:np.pi:complex(Num)]

        Para, Perp = GetFieldsFromMesh(m                    = Parent.Index,
                                       x                    = Parent.SizeParam,
                                       ThetaMesh            = self['Theta'].flatten(),
                                       PhiMesh              = self['Phi'].flatten()-np.pi/2,
                                       Polarization         = 0)

        self['Parallel'], self['Perpendicular'] = Para, Perp
        self['SPF'] = np.sqrt( Para.__abs__()**2 + Perp.__abs__()**2 )


    def Plot(self):

        x, y, z = cs.sp2cart(self['SPF'].reshape(self['Phi'].shape), self['Phi'], self['Theta'])

        fig, ax = plt.subplots(1, figsize=(3, 3), subplot_kw = {'projection':'3d'})
        ax.set_title(r'Complex Scattering Phase Function: Real{$ E_{||}$}')
        ax.set_ylabel(r'Y-direction')
        ax.set_xlabel(r'X-direction')
        ax.set_zlabel(r'Z-direction')

        scamap = plt.cm.ScalarMappable(cmap='jet')

        ax.plot_surface( x, y, z,
                         rstride     = 2,
                         cstride     = 2,
                         linewidth   = 0.3,
                         facecolors  = scamap.to_rgba(self['Parallel'].real.reshape(self['Phi'].shape) ),
                         alpha =1,
                         edgecolors  = 'k',
                         #antialiased = False,
                         #shade  =False
                         )

        xLim = ax.get_xlim(); yLim = ax.get_ylim(); zLim = ax.get_zlim()

        Min = min(xLim[0],yLim[0]); Max = max(xLim[1], yLim[1])

        ax.set_xlim(Min, Max); ax.set_ylim(Min, Max)

        plt.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            SPF, Parallel, Perpendicular, Theta, Phi
        Structured data: Yes
        Shape:           {self['Phi'].shape}"""



class S1S2(dict):
    """Short summary.

    Parameters
    ----------
    SizeParam : np.array
        Description of parameter `SizeParam`.
    Index : float
        Description of parameter `Index`.
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

    def __init__(self, Parent, Num):
        self['Phi'] = np.linspace(0, 2*np.pi, Num)
        self['S1'], self['S2'] = GetS1S2(Parent.Index, Parent.SizeParam, self['Phi'])


    def Plot(self) -> None:


        S1 = np.abs(self['S1'])
        S2 = np.abs(self['S2'])

        fig, axes = plt.subplots(nrows      = 1,
                                 ncols      = 2,
                                 subplot_kw = {'projection':'polar'})

        axes[0].set_title('S1 function'); axes[1].set_title('S2 function')

        axes[0].plot(self['Phi'],
                     S1,
                     color = 'k')

        axes[0].fill_between(x     = self['Phi'],
                             y2    = 0,
                             y1    = S1,
                             color = 'C0',
                             alpha = 0.4)


        axes[1].plot(self['Phi'],
                     S2,
                     color = 'k')

        axes[1].fill_between(x     = self['Phi'],
                             y2    = 0,
                             y1    = S2,
                             color = 'C1',
                             alpha = 0.4)

        plt.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            S1, S2, Phi
        Structured data: Yes
        Shape:           {self['Phi'].shape}"""




class ScalarFarField(dict):


    def __init__(self, Num = 200, Parent = None):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(Num), -np.pi/2:np.pi/2:complex(Num)]

        Para, Perp = GetFieldsFromMesh(m            = Parent.Index,
                                       x            = Parent.SizeParam,
                                       ThetaMesh    = Theta.flatten(),
                                       PhiMesh      = Phi.flatten() - np.pi/2,
                                       Polarization = Parent.Source.Polarization.Radian)

        self['Parallel'] = Para

        self['Perpendicular'] = Perp

        self['Theta'] = Theta

        self['Phi'] = Phi


    def Plot(self, num=200):

        PlotStructuredSphere(Phi     = np.rad2deg(self['Phi']),
                             Theta   = np.rad2deg(self['Theta']),
                             Scalar  = self['Parallel'].reshape(self['Theta'].shape))

        PlotStructuredSphere(Phi     = np.rad2deg(self['Phi']),
                             Theta   = np.rad2deg(self['Theta']),
                             Scalar  = self['Perpendicular'].reshape(self['Theta'].shape))

    def __repr__(self):

        return f"""
        Object:          Dictionary
        Keys:            Parallel, Perpendicular, Theta, Phi
        Structured data: Yes
        Shape:           {self['Theta'].shape}"""







# -
