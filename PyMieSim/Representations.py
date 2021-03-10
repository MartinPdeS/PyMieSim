import numpy as np
from ai import cs
from mayavi import mlab
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"


from PyMieSim.utils import PlotStructuredSphere, LoadLibraries
from PyMieSim.utils import Direct2spherical, AngleUnit2DirectUnit
from PyMieSim.LMT.Sphere import S1S2 as LMTS1S2, FieldsUnstructured as LMTFields
from PyMieSim.GLMT.Sphere import S1S2Structured as GLMTS1S2, FieldsStructured as GLMTFieldsStructured

class Stokes(np.ndarray):
    """ http://adsabs.harvard.edu/pdf/1983AuJPh..36..701B

    """

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

        axes[0].pcolormesh(self.Parent.Meshes.Theta.Radian,
                           -( self.Parent.Meshes.Phi.Radian - np.pi/2),
                           np.real(self.Perpendicular),
                           shading='auto')

        axes[1].pcolormesh(self.Parent.Meshes.Theta.Radian,
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

        phi, theta = np.linspace(-np.pi/2, np.pi/2, Num), np.linspace(-np.pi, np.pi, Num)
        self['Phi'], self['Theta'] = np.meshgrid(phi, theta)


        if Parent.Source.GLMT:
            Para, Perp = GLMTFieldsStructured(Index        = Parent.Index,
                                              Diameter     = Parent.Diameter,
                                              Wavelength   = Parent.Source.Wavelength,
                                              nMedium      = Parent.nMedium,
                                              Phi          = phi-np.pi/2,
                                              Theta        = theta,
                                              Polarization = Parent.Source.Polarization.Radian,
                                              E0           = float(Parent.Source.E0),
                                              R            = 1.,
                                              BSC          = Parent.Source._BSC_,
                                              MaxOrder     = Parent.Source.MaxOrder)

            self['EPhi'], self['ETheta'] = Para, Perp

            self['SPF'] = np.sqrt( Para.__abs__()**2 + Perp.__abs__()**2 ).reshape([Num, Num]).T


        else:
            Para, Perp = LMTFields(Index        = Parent.Index,
                                   Diameter     = Parent.Diameter,
                                   Wavelength   = Parent.Source.Wavelength,
                                   nMedium      = 1.0,
                                   Phi          = self['Phi'].flatten()-np.pi/2,
                                   Theta        = self['Theta'].flatten(),
                                   Polarization = 0.,
                                   E0           = Parent.Source.E0,
                                   R            = 1.,
                                   Lenght       = self['Theta'].flatten().size)

            self['EPhi'], self['ETheta'] = Para, Perp

            self['SPF'] = np.sqrt( Para.__abs__()**2 + Perp.__abs__()**2 ).reshape([Num, Num])


    def Plot(self):

        x, y, z = cs.sp2cart(self['SPF'], self['Phi'], self['Theta'])

        X = np.linspace(np.min(x),np.max(x),10)*1.3
        Y = np.linspace(np.min(y),np.max(y),10)*1.3
        Z = np.linspace(np.min(z),np.max(z),10)*1.3

        radius = np.abs( min(np.min(x), np.min(y), np.min(z))/100 )

        mlab.plot3d(X,
                    np.zeros_like(X),
                    np.zeros_like(X),
                    line_width=1e-12,
                    tube_radius=radius
                    )

        mlab.plot3d(np.zeros_like(Y),
                    Y,
                    np.zeros_like(Y),
                    line_width=1e-12,
                    tube_radius=radius
                    )

        mlab.plot3d(np.zeros_like(Z),
                    np.zeros_like(Z),
                    Z,
                    line_width=1e-12,
                    tube_radius=radius
                    )

        mlab.mesh(x, y, z, colormap='viridis')

        mlab.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            SPF, EPhi, ETheta, Theta, Phi
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""



class S1S2(dict):

    def __init__(self, Parent, Num):

        self['Phi'] = np.linspace(0, 2*np.pi, Num)

        Phi = np.linspace(-np.pi,np.pi,100);

        args = (1.8, 3e-6, 1e-6, 1, Phi)
        self['S1'], self['S2'] = LMTS1S2(Index      = Parent.Index,
                                         Diameter   = Parent.Diameter,
                                         Wavelength = Parent.Source.Wavelength,
                                         nMedium    = Parent.nMedium,
                                         Phi        = self['Phi'])


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
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""




class ScalarFarField(dict):

    def __init__(self, Num = 200, Parent = None, Distance=1):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(Num), -np.pi/2:np.pi/2:complex(Num)]

        EPhi, ETheta = LMTFields(Index        = Parent.Index,
                                 Diameter     = Parent.Diameter,
                                 Wavelength   = Parent.Source.Wavelength,
                                 nMedium      = Parent.nMedium,
                                 Phi          = Phi.flatten() - np.pi/2,
                                 Theta        = Theta.flatten(),
                                 Polarization = Parent.Source.Polarization.Radian,
                                 E0           = Parent.Source.E0,
                                 R            = Distance,
                                 Lenght       = Phi.flatten().size)


        self['Distance'] = Distance

        self['Theta'] = Theta

        self['Phi'] = Phi

        self['EPhi'] = EPhi.reshape(self['Theta'].shape)

        self['ETheta'] = ETheta.reshape(self['Theta'].shape)



    def Plot(self):
        """Method plots the scattered Far-Field
        :math:`E_{\\phi}(\\phi,\\theta)^2 , E_{\\theta}(\\phi,\\theta)^2`.

        Parameters
        ----------
        Num : :class:`int`
            Number of point to spatially (:math:`\\theta , \\phi`) evaluate the SPF [Num, Num].

        """

        fig, axes = plt.subplots(2,2,figsize=(6,4),subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea()})

        for F, Name in enumerate( ['EPhi', 'ETheta'] ):

            for m, Part in enumerate(['Real', 'Imaginary']) :

                if m == 0: data = self[Name].real
                if m == 1: data = self[Name].imag

                im0 = axes[m,F].contourf(np.rad2deg(self['Theta']),
                                         np.rad2deg(self['Phi']),
                                         data,
                                         antialiased=False,
                                         cmap = 'inferno',
                                         #levels = 100,
                                         transform = ccrs.PlateCarree())

                gl = axes[m,F].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
                gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
                gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

                plt.colorbar(mappable=im0, fraction=0.046, orientation='vertical', ax=axes[m,F])

                if F == 0: name = r"$E_{\phi}(\phi, \theta)$"
                if F == 1: name = r"$E_{\theta}(\phi, \theta)$"

                axes[m,F].set_title(f'{Part} Part' + '\n' +  name, fontsize=8)

        fig.tight_layout()

        plt.show()



    def __repr__(self):

        return f"""
        Object:          Dictionary
        Keys:            EPhi, ETheta, Theta, Phi, Distance
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Theta'].shape}"""





class Footprint(dict):

    def __init__(self, Scatterer, Detector, Num=100):

        x, y = np.mgrid[-50: 50: complex(Num), -50: 50: complex(Num)]

        MaxAngle = np.abs( np.pi/2 - Detector.Mesh.Phi.Radian.min() )

        Phi, Theta = Direct2spherical(X=x, Y=y, MaxAngle=MaxAngle)

        Direct = AngleUnit2DirectUnit(Phi, Scatterer.Source.k)

        Perp =  (Detector.StructuredFarField(Num=Num, SFactor=16) *\
        Scatterer.Perpendicular( Phi.flatten(), Theta.flatten() ).reshape(Theta.shape) )

        Para = (Detector.StructuredFarField(Num=Num, SFactor=16) *\
        Scatterer.Parallel( Phi.flatten(), Theta.flatten() ).reshape(Theta.shape) )

        n = 5

        FourierPara = np.fft.ifft2(Para, s=[512*n, 512*n])

        FourierPara = np.fft.fftshift(FourierPara).__abs__()**2

        FourierPerp = np.fft.ifft2(Perp,  s=[512*n, 512*n])

        FourierPerp = np.fft.fftshift(FourierPerp).__abs__()**2

        Direct = np.linspace(np.min(Direct), np.max(Direct), FourierPerp.shape[0])

        self['Map'] = FourierPara + FourierPerp

        self['DirectX'] = Direct; self['DirectY'] = Direct



    def Plot(self):

        fig = plt.figure()

        ax = fig.add_subplot(111)

        ax.contourf(self['DirectX']*1e6, self['DirectY']*1e6, self['Map'], cmap='gray')

        ax.set_xlabel(r'Offset distance in X-axis [$\mu$m]')

        ax.set_ylabel(r'Offset distance in Y-axis [$\mu$m]')

        ax.set_title('Scatterer Footprint')

        plt.show()


    def __repr__(self):
        return f"""
        Object:          Dictionary
        Keys:            Map, DirectX, DirectY
        Structured data: Yes
        Method:          <Plot>
        Shape:           {self['Phi'].shape}"""



# -
