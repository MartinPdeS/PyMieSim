import numpy as np
from ai import cs
from mayavi import mlab
from scipy.interpolate import griddata

from PyMieSim.Physics import Angle
import PyMieSim


def GetFieldBinding(Scatterer, Structured, R, Phi, Theta):

    kwarg = { 'Index'       : Scatterer.Index,
               'Diameter'   : Scatterer.Diameter,
               'Wavelength' : Scatterer.Source.Wavelength,
               'nMedium'    : Scatterer.nMedium,
               'Phi'        : Phi,
               'Theta'      : Theta,
               'R'          : R,
               'E0'         : Scatterer.Source.E0}


    lib = PyMieSim
    if Scatterer.Source.GLMT:
        lib = lib.GLMT
        kwarg['BSC'] = Scatterer.Source._BSC_
        kwarg['MaxOrder'] = Scatterer.Source.MaxOrder

    else:
        lib = lib.LMT

    if Scatterer.type == 'Sphere': lib = lib.Sphere

    if Scatterer.type == 'Cylinder': lib = lib.Cylinder

    if Structured:
        lib = lib.Structured
    else:
        lib = lib.Unstructured

    if Scatterer.Source.Polarization:
        return lib.Fields(**kwarg, Polarization = Scatterer.Source.Polarization.Radian )
    else:
        return lib.FieldsUnpolarized(**kwarg)



def _GetFieldBinding(Scatterer, Structured, R, Phi, Theta):

    kwarg = { 'Index'       : Scatterer.Index,
               'Diameter'   : Scatterer.Diameter,
               'Wavelength' : Scatterer.Source.Wavelength,
               'nMedium'    : Scatterer.nMedium,
               'Phi'        : Phi,
               'Theta'      : Theta,
               'R'          : R,
               'E0'         : Scatterer.Source.E0}

    if Structured:
        if Scatterer.Source.GLMT:
            if Scatterer.Source.Polarization:
                from PyMieSim.GLMT.Sphere import FieldsStructured
                return FieldsStructured(**kwarg,
                                          Polarization = Scatterer.Source.Polarization.Radian,
                                          BSC          = Scatterer.Source._BSC_,
                                          MaxOrder     = Scatterer.Source.MaxOrder)

            else:
                from PyMieSim.GLMT.Sphere import FieldsStructuredUnpolarized
                return FieldsStructuredUnpolarized(**kwarg,
                                                     BSC          = Scatterer.Source._BSC_,
                                                     MaxOrder     = Scatterer.Source.MaxOrder)


        else:
            if Scatterer.Source.Polarization:
                from PyMieSim.LMT.Sphere import FieldsStructured
                return FieldsStructured(**kwarg,
                                          Polarization = Scatterer.Source.Polarization.Radian)
            else:
                from PyMieSim.LMT.Sphere import FieldsStructuredUnpolarized
                return FieldsStructuredUnpolarized(**kwarg)

    else:

        if Scatterer.Source.GLMT:
            if Scatterer.Source.Polarization:
                from PyMieSim.GLMT.Sphere import FieldsUnstructured
                return FieldsUnstructured(**kwarg,
                                            Polarization = Scatterer.Source.Polarization.Radian,
                                            BSC          = Scatterer.Source._BSC_,
                                            MaxOrder     = Scatterer.Source.MaxOrder)

            else:
                from PyMieSim.GLMT.Sphere import FieldsUnstructuredUnpolarized
                return FieldsUnstructuredUnpolarized(**kwarg,
                                                       BSC          = Scatterer.Source._BSC_,
                                                       MaxOrder     = Scatterer.Source.MaxOrder)


        else:
            if Scatterer.Source.Polarization:
                from PyMieSim.LMT.Sphere import FieldsUnstructured
                return FieldsUnstructured(**kwarg,
                                            Polarization = Scatterer.Source.Polarization.Radian)
            else:
                from PyMieSim.LMT.Sphere import FieldsUnstructuredUnpolarized
                return FieldsUnstructuredUnpolarized(**kwarg)




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

    grid = griddata((x, y), v, (xp, yp), method=algorithm).ravel()

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
    v[nans] = griddata((x[notnans], y[notnans]),
                        v[notnans],
                        (x[nans], y[nans]),
                        method='nearest').ravel()
    return v





def UnitSphere(Num, Radius=1.):
    pi = np.pi
    phi, theta = np.mgrid[-pi/2:pi/2:complex(Num), -pi:pi:complex(Num)]
    x, y, z = cs.sp2cart(phi*0+Radius, phi, theta)

    return x, y, z


def PlotUnstructured(Scalar, Mesh, Name=''):
    x, y, z = cs.sp2cart(Mesh.Phi.Radian.flatten()*0+1,
                         Mesh.Phi.Radian.flatten(),
                         Mesh.Theta.Radian.flatten())

    offset = 3
    X = np.linspace(-1,1,10)*1.3
    Y = np.linspace(-1,1,10)*1.3
    Z = np.linspace(-1,1,10)*2

    zeros = np.zeros_like(X)

    fig = mlab.figure(size=(600,300))

    mlab.plot3d(X, zeros, zeros, line_width=1e-12)
    mlab.plot3d(zeros, Y, zeros, line_width=1e-12)
    mlab.plot3d(zeros, zeros, Z, line_width=1e-12)

    mlab.plot3d(X+offset, zeros, zeros, line_width=1e-12)
    mlab.plot3d(zeros+offset, Y, zeros, line_width=1e-12)
    mlab.plot3d(zeros+offset, zeros, Z, line_width=1e-12)

    mlab.text(0, 0, u'Z', z=2, width=0.01)
    mlab.text(0, 1.3, u'Y', z=0, width=0.01)
    mlab.text(1.3, 0, u'X', z=0, width=0.01)

    mlab.text(0+offset, 0, u'Z', z=2, width=0.01)
    mlab.text(0+offset, 1.3, u'Y', z=0, width=0.01)
    mlab.text(1.3+offset, 0, u'X', z=0, width=0.01)

    mlab.text(offset, 0, u'Imaginary', z=-2., width=0.15)
    mlab.text(0, 0, u'Real', z=-2, width=0.07)

    mlab.text(offset/2, 0, Name, z=2, width=0.2)

    xp, yp, zp = UnitSphere(Num=50, Radius=1.)

    mlab.mesh(xp, yp, zp, colormap='gray', opacity=0.5)

    mlab.mesh(xp+offset, yp, zp, colormap='gray', opacity=0.5)

    im0 = mlab.points3d(x, y, z, Scalar.real, mode='sphere', scale_mode='none', colormap='inferno')

    im1 = mlab.points3d(x+offset, y, z, Scalar.imag, mode='sphere', scale_mode='none', colormap='inferno')

    mlab.colorbar(object = im0, label_fmt="%.0e", nb_labels=5, title='Real part', orientation='horizontal' )

    mlab.colorbar(object = im1, label_fmt="%.0e", nb_labels=5, title='Imaginary part', orientation='vertical' )

    mlab.show()







def PlotStructuredAbs(Scalar, Phi, Theta, Name=''):

        x, y, z = cs.sp2cart(Scalar, Phi, Theta)

        X = np.linspace(np.min(x),np.max(x),10)*1.3
        Y = np.linspace(np.min(y),np.max(y),10)*1.3
        Z = np.linspace(np.min(z),np.max(z),10)*1.3

        radius = np.abs( min(np.min(x), np.min(y), np.min(z))/100 )

        zeros = np.zeros_like(X)

        mlab.plot3d(X, zeros, zeros, line_width=1e-12, tube_radius=radius)
        mlab.plot3d(zeros, Y, zeros, line_width=1e-12, tube_radius=radius)
        mlab.plot3d(zeros, zeros, Z, line_width=1e-12, tube_radius=radius)

        mlab.text(0, 0, u'Z', z=Z[-1], width=0.01)
        mlab.text(0, Y[-1], u'Y', z=0, width=0.01)
        mlab.text(X[-1], 0, u'X', z=0, width=0.01)

        mlab.text(0, 0, Name, z=Z[-1]*1.15, width=0.5)

        mlab.mesh(x, y, z, colormap='viridis')

        mlab.show()


def PlotStructuredAmplitude(Scalar, Phi, Theta, Name=''):
    x, y, z = cs.sp2cart(Phi*0+1,
                         Phi,
                         Theta)

    offset = 3
    X = np.linspace(-1,1,10)*1.3
    Y = np.linspace(-1,1,10)*1.3
    Z = np.linspace(-1,1,10)*2

    zeros = np.zeros_like(X)

    fig = mlab.figure(size=(600,300))

    mlab.plot3d(X, zeros, zeros, line_width=1e-12)
    mlab.plot3d(zeros, Y, zeros, line_width=1e-12)
    mlab.plot3d(zeros, zeros, Z, line_width=1e-12)

    mlab.plot3d(X+offset, zeros, zeros, line_width=1e-12)
    mlab.plot3d(zeros+offset, Y, zeros, line_width=1e-12)
    mlab.plot3d(zeros+offset, zeros, Z, line_width=1e-12)

    mlab.text(0, 0, u'Z', z=2, width=0.01)
    mlab.text(0, 1.3, u'Y', z=0, width=0.01)
    mlab.text(1.3, 0, u'X', z=0, width=0.01)

    mlab.text(0+offset, 0, u'Z', z=2, width=0.01)
    mlab.text(0+offset, 1.3, u'Y', z=0, width=0.01)
    mlab.text(1.3+offset, 0, u'X', z=0, width=0.01)

    mlab.text(offset, 0, u'Imaginary', z=-2., width=0.15)
    mlab.text(0, 0, u'Real', z=-2, width=0.07)

    mlab.text(offset/2, 0, Name, z=2, width=0.2)

    xp, yp, zp = UnitSphere(Num=50, Radius=1.)

    mlab.mesh(xp, yp, zp, colormap='gray', opacity=0.5)

    mlab.mesh(xp+offset, yp, zp, colormap='gray', opacity=0.5)

    im0 = mlab.points3d(x, y, z, Scalar.real, mode='sphere', scale_mode='none', colormap='inferno')

    im1 = mlab.points3d(x+offset, y, z, Scalar.imag, mode='sphere', scale_mode='none', colormap='inferno')

    mlab.colorbar(object = im0, label_fmt="%.0e", nb_labels=5, title='Real part', orientation='horizontal' )

    mlab.colorbar(object = im1, label_fmt="%.0e", nb_labels=5, title='Imaginary part', orientation='vertical' )

    mlab.show()


""" Ref: https://optiwave.com/optifdtd-manuals/fdtd-far-field-transform/"""




def Angle2Direct(AngleVec: np.ndarray, k: float,) -> np.ndarray:

    RadSpace = np.deg2rad(AngleVec)

    FourierSpace = np.sin(RadSpace) * k / (2 * np.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = np.fft.fftshift( np.fft.fftfreq( AngleVec.shape[0], d = fourier_unit ) )

    return DirectSpace



def Direct2Angle(DirectVec: np.ndarray, k: float) -> np.ndarray:

    direct_unit = (DirectVec[1] - DirectVec[0]).__abs__()

    FourierSpace = np.fft.fftshift( np.fft.fftfreq( DirectVec.shape[0], d = direct_unit ) )

    AngleVec = np.arcsin(2 * np.pi * FourierSpace / k) # conversion spatial frequency to angular space

    if np.isnan(AngleVec).any():
        raise Exception("Magnification too large.")

    return AngleVec * 180 / np.pi



def NA2Angle(NA: float) -> np.ndarray:

    return Angle( np.arcsin(NA), unit='Radian')



def Direct2spherical(X, Y, MaxAngle):
    Z = 50 / np.tan(MaxAngle)

    _, Phi, Theta = cs.cart2sp(X, Y, X*0+Z)

    return Phi, Theta

def Direct2Angle(X, Y, MaxAngle):
    MaxZ = np.max(X) / np.cos(MaxAngle)



def AngleUnit2DirectUnit(Angle, k):
    FourierSpace = np.sin(Angle) * k / (2 * np.pi)

    fourier_unit = (FourierSpace[1] - FourierSpace[0]).__abs__()

    DirectSpace = np.fft.fftshift( np.fft.fftfreq( Angle.shape[0], d = fourier_unit ) )

    return DirectSpace
# -
