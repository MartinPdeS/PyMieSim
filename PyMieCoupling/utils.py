import numpy as np
import fibermodes
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import scipy


from PyMieCoupling.cpp.S1S2 import GetFieldsFromMesh

class Source(object):

    def __init__(self,
                 Wavelength:   float,
                 Polarization: float,
                 Power:        float = 1,
                 Radius:       float = 1):

        self.Wavelength = Wavelength

        self.k = 2 * np.pi / Wavelength

        self.Power = Power

        if Polarization != None:
            self.Polarization = _Polarization(Polarization)
        else:
            self.Polarization = None


class _Polarization(object):

    def __init__(self, input,):
        if input == 'None':
            self.Degree = 'None'
            self.Radian = 'None'
        else:
            self.Degree = input
            self.Radian = np.deg2rad(input)



class Angle(object):

    def __init__(self, input, unit='Degree'):
        if input == 'None':
            self.Degree = 'None'
            self.Radian = 'None'

        if unit == 'Degree':
            self.Degree = input
            self.Radian = np.deg2rad(input)
        if unit == 'Radian':
            self.Degree = np.rad2deg(input)
            self.Radian = input


def SMF28():
    CoreDiameter = 8.2e-6
    cladDiameter = 125e-6

    Fiber = fiber(core_radius = CoreDiameter,
                  core_index  = 1.4456,
                  clad_radius = cladDiameter,
                  clad_index  = 1.4444)

    return Fiber, CoreDiameter



class fiber(object):

    def __init__(self,
                 core_radius,
                 core_index,
                 clad_radius,
                 clad_index):

        self.MaxDirect = 2 * clad_radius

        factory = fibermodes.FiberFactory()

        factory.addLayer(name     = 'core',
                         radius   = core_radius,
                         material = 'Fixed',
                         geometry = "StepIndex",
                         index    = 1.4489)

        factory.addLayer(name     = 'cladding',
                         material = 'Fixed',
                         index    = 1)

        self.source = factory[0]



def InterpFull(Meshes, Scalar, Shape):

    Phi, Theta = np.mgrid[-np.pi/2:np.pi/2:complex(Shape[0]),
                          -np.pi:np.pi:complex(Shape[1])]

    Scalar = interp_at(Meshes.Phi.Radian,
                       Meshes.Theta.Radian,
                       Scalar.astype(np.complex).flatten(),
                       Phi.flatten(),
                       Theta.flatten(),
                       algorithm='linear',
                       extrapolate=True)

    return Scalar.reshape(Shape), Phi, Theta



def PlotUnstructureData(Scalar, phi, theta):

    fig, ax = plt.subplots(1,2,figsize=(15,8))

    im0 = ax[0].tripcolor(theta, phi, Scalar.real)

    im1 = ax[1].tripcolor(theta, phi, Scalar.imag)

    plt.colorbar(mappable=im0, ax=ax[0])

    plt.colorbar(mappable=im1, ax=ax[1])

    ax[0].plot(theta, phi, 'ko ', markersize=2)

    ax[1].plot(theta, phi, 'ko ', markersize=2)

    plt.show()



def interp_at(x, y, v, xp, yp, algorithm='cubic', extrapolate=False):
    """
    Interpolate data onto the specified points.

    Parameters:

    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.
    * xp, yp : 1D arrays
        Points where the data values will be interpolated
    * algorithm : string
        Interpolation algorithm. Either ``'cubic'``, ``'nearest'``,
        ``'linear'`` (see scipy.interpolate.griddata)
    * extrapolate : True or False
        If True, will extrapolate values outside of the convex hull of the data
        points.

    Returns:

    * v : 1D array
        1D array with the interpolated v values.

    """
    if algorithm not in ['cubic', 'linear', 'nearest']:
        raise ValueError("Invalid interpolation algorithm: " + str(algorithm))

    grid = scipy.interpolate.griddata((x, y),
                                      v,
                                      (xp, yp),
                                      method=algorithm).ravel()

    if extrapolate and algorithm != 'nearest' and np.any(np.isnan(grid)):
        grid = extrapolate_nans(xp, yp, grid)
    return grid


def extrapolate_nans(x, y, v):
    """
    Extrapolate the NaNs or masked values in a grid INPLACE using nearest
    value.

    .. warning:: Replaces the NaN or masked values of the original array!

    Parameters:

    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.

    Returns:

    * v : 1D array
        The array with NaNs or masked values extrapolated.

    """
    if np.ma.is_masked(v):
        nans = v.mask
    else:
        nans = np.isnan(v)
    notnans = np.logical_not(nans)
    v[nans] = scipy.interpolate.griddata((x[notnans], y[notnans]), v[notnans],
                                         (x[nans], y[nans]),
                                         method='nearest').ravel()
    return v




def PlotFarField(Phi, Theta, Scalar, Meshes=None, Name='Field', scatter=True):

        fig, axes = plt.subplots(nrows      = 1,
                                  ncols      = 2,
                                  figsize    = (8, 4),
                                  subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea()})

        im0 = axes[0].contourf(np.rad2deg(Theta),
                               np.rad2deg(Phi),
                               Scalar.real,
                               cmap      = 'inferno',
                               transform = ccrs.PlateCarree(),
                               levels    = 20
                               )

        gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
        gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

        plt.colorbar(mappable=im0, fraction=0.046, orientation='vertical', ax=axes[0])
        axes[0].set_title('Real Part {}'.format(Name))
        axes[0].set_ylabel(r'Angle $\phi$ [Degree]')
        axes[0].set_xlabel(r'Angle $\theta$ [Degree]')
        axes[0].grid(True, which='minor')


        im1 = axes[1].contourf(np.rad2deg(Theta),
                               np.rad2deg(Phi),
                               Scalar.imag,
                               cmap      = 'inferno',
                               transform = ccrs.PlateCarree(),
                               levels    = 20
                               )

        gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
        gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
        gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

        plt.colorbar(mappable=im1, fraction=0.046, orientation='vertical', ax=axes[1])
        axes[1].set_title('Imaginary Part {}'.format(Name))
        axes[1].set_ylabel(r'Angle $\phi$ [Degree]')
        axes[1].set_xlabel(r'Angle $\theta$ [Degree]')
        axes[1].grid(True, which='minor')


        fig.tight_layout()

        plt.show()


        return fig



def PlotStructuredSphere(Scalar, Phi, Theta, Name=''):

    fig, axes = plt.subplots(1,2,figsize=(8,4),subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea()})

    im0 = axes[0].contourf(Theta,
                            Phi,
                            Scalar.real,
                            antialiased=False,
                            cmap = 'inferno',
                            #levels = 100,
                            transform = ccrs.PlateCarree())

    im1 = axes[1].contourf(Theta,
                            Phi,
                            Scalar.imag,
                            antialiased=False,
                            cmap = 'inferno',
                            #levels = 100,
                            transform = ccrs.PlateCarree())



    gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

    gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

    plt.colorbar(mappable=im0, fraction=0.046, orientation='vertical', ax=axes[0])
    plt.colorbar(mappable=im1, fraction=0.046, orientation='vertical', ax=axes[1])

    axes[0].set_title(f'Real Part {Name}')
    axes[0].set_ylabel(r'Angle $\phi$ [Degree]')
    axes[0].set_xlabel(r'Angle $\theta$ [Degree]')


    axes[1].set_title(f'Imaginary Part {Name}')
    axes[1].set_ylabel(r'Angle $\phi$ [Degree]')
    axes[1].set_xlabel(r'Angle $\theta$ [Degree]')


    #axes[0].set_extent([-12755636.1863, 12755636.1863, -12727770.598700099, 12727770.598700099])
    #axes[1].set_extent([-170, 170, -90, 90], crs=ccrs.PlateCarree())

    fig.tight_layout()

    plt.show()



def PlotUnstructuredSphere(Scalar, Phi, Theta, Name=''):

    fig, axes = plt.subplots(1,
                             2,
                             figsize=(8,4),
                             subplot_kw = {'projection':ccrs.LambertAzimuthalEqualArea()}
    )

    im0 = axes[0].tripcolor(Theta,
                            Phi,
                            Scalar.real,
                            antialiased=False,
                            transform = ccrs.PlateCarree()
                            )

    im1 = axes[1].tricontour(Theta,
                              Phi,
                              Scalar.imag,
                              antialiased=False,
                              transform = ccrs.PlateCarree()
                             )



    gl = axes[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

    gl = axes[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(-180,181,30))
    gl.ylocator = matplotlib.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

    plt.colorbar(mappable=im0, fraction=0.046, orientation='vertical', ax=axes[0])
    plt.colorbar(mappable=im1, fraction=0.046, orientation='vertical', ax=axes[1])

    axes[0].set_title(f'Real Part {Name}')
    axes[0].set_ylabel(r'Angle $\phi$ [Degree]')
    axes[0].set_xlabel(r'Angle $\theta$ [Degree]')
    axes[0].grid(True, which='minor')

    axes[1].set_title(f'Imaginary Part {Name}')
    axes[1].set_ylabel(r'Angle $\phi$ [Degree]')
    axes[1].set_xlabel(r'Angle $\theta$ [Degree]')
    axes[1].grid(True, which='minor')


    axes[0].set_extent([-170, 170, -90, 90], crs=ccrs.PlateCarree())
    axes[1].set_extent([-170, 170, -90, 90], crs=ccrs.PlateCarree())


    axes[0].plot(Theta,
                 Phi,
                 'ko ',
                 markersize=0.2,
                 transform = ccrs.PlateCarree()
                 )

    axes[1].plot(Theta,
                 Phi,
                 'ko ',
                 markersize=0.2,
                 transform = ccrs.PlateCarree()
                 )
    fig.tight_layout()


    plt.show()
