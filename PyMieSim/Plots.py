import numpy as np
pi = np.pi
from mayavi import mlab
from scipy.interpolate import griddata

from PyMieSim.Physics import Angle
from PyMieSim.utils import Sp2Cart, Cart2Sp



def UnitAxes(Figure, Radius=1., xOffset=0.):
    X = np.linspace(-1.0,1.0,10)*1.3 * Radius
    Y = np.linspace(-1.0,1.0,10)*1.3 * Radius
    Z = np.linspace(-1.0,1.0,10)*2.0 * Radius

    zeros = np.zeros_like(X)

    mlab.plot3d(X + xOffset, zeros, zeros, line_width=1e-12, figure=Figure)
    mlab.plot3d(zeros + xOffset, Y, zeros, line_width=1e-12, figure=Figure)
    mlab.plot3d(zeros + xOffset, zeros, Z, line_width=1e-12, figure=Figure)

    mlab.text(0.0 + xOffset, 0.0, u'Z', z=2.0, width=0.01, figure=Figure)
    mlab.text(0.0 + xOffset, 1.3, u'Y', z=0.0, width=0.01, figure=Figure)
    mlab.text(1.3 + xOffset, 0.0, u'X', z=0.0, width=0.01, figure=Figure)


def UnitSphere(Num, Radius=1.):

    phi, theta = np.mgrid[-pi/2:pi/2:complex(Num), -pi:pi:complex(Num)]
    x, y, z = Sp2Cart(phi*0+Radius, phi, theta)

    return x, y, z



def PlotUnstructuredAbs(Scalar, Mesh, Name=''):

    x, y, z = Sp2Cart(Mesh.Phi.Radian.flatten()*0+1,
                      Mesh.Phi.Radian.flatten(),
                      Mesh.Theta.Radian.flatten())

    _PlotUnstructuredAbs(Scalar  = Mesh.Phi.Radian*0+1,
                         Phi     = Mesh.Phi.Radian,
                         Theta   = Mesh.Theta.Radian,
                         Name    = Name)




def _PlotUnstructuredAbs(Scalar, Phi, Theta, Name=''):

    Phi = Phi.flatten()
    Theta = Theta.flatten()
    Scalar = Scalar.flatten()

    x, y, z = Sp2Cart(Phi*0+1,
                      Phi,
                      Theta)

    offset = 3
    X = np.linspace(-1,1,10)*1.3
    Y = np.linspace(-1,1,10)*1.3
    Z = np.linspace(-1,1,10)*2

    zeros = np.zeros_like(X)

    fig = mlab.figure(size=(600,300))

    mlab.plot3d(X+offset, zeros, zeros, line_width=1e-12)
    mlab.plot3d(zeros+offset, Y, zeros, line_width=1e-12)
    mlab.plot3d(zeros+offset, zeros, Z, line_width=1e-12)


    mlab.text(0.0 + offset, 0.0, u'Z', z=2.0, width=0.01)
    mlab.text(0.0 + offset, 1.3, u'Y', z=0.0, width=0.01)
    mlab.text(1.3 + offset, 0.0, u'X', z=0.0, width=0.01)

    mlab.text(0 + offset, 0, Name, z=2, width=0.2)

    xp, yp, zp = UnitSphere(Num=50, Radius=1.)

    mlab.mesh(xp + offset, yp, zp, colormap='gray', opacity=0.5)

    im0 = mlab.points3d(x + offset, y, z, np.abs(Scalar), mode='sphere', scale_mode='none', colormap='inferno')

    mlab.colorbar(object = im0, label_fmt="%.0e", nb_labels=5, title='Real part', orientation='horizontal' )



def PlotUnstructured(Scalar, Mesh, Name=''):
    x, y, z = Sp2Cart(Mesh.Phi.Radian.flatten()*0+1,
                      Mesh.Phi.Radian.flatten(),
                      Mesh.Theta.Radian.flatten())

    offset = 3
    X = np.linspace(-1,1,10)*1.3
    Y = np.linspace(-1,1,10)*1.3
    Z = np.linspace(-1,1,10)*2

    zeros = np.zeros_like(X)

    fig = mlab.figure(size=(600,300))

    UnitAxes(Figure=fig, Radius=1., xOffset=0.)

    UnitAxes(Figure=fig, Radius=1., xOffset=offset)

    mlab.text(offset, 0, u'Imaginary', z=-2., width=0.15)
    mlab.text(0, 0, u'Real', z=-2, width=0.07)

    mlab.text(offset/2, 0, Name, z=2, width=0.2)

    xp, yp, zp = UnitSphere(Num=50, Radius=1.)

    mlab.mesh(xp, yp, zp, colormap='gray', opacity=0.5)

    mlab.mesh(xp+offset, yp, zp, colormap='gray', opacity=0.5)

    im0 = mlab.points3d(*Mesh.CartCoord, Scalar.real, mode='sphere', scale_mode='none', colormap='inferno')

    im1 = mlab.points3d(Mesh.CartCoord[0]+offset, Mesh.CartCoord[1], Mesh.CartCoord[2], Scalar.imag, mode='sphere', scale_mode='none', colormap='inferno')

    mlab.colorbar(object = im0, label_fmt="%.0e", nb_labels=5, title='Real part', orientation='horizontal' )

    mlab.colorbar(object = im1, label_fmt="%.0e", nb_labels=5, title='Imaginary part', orientation='vertical' )


def PlotStructuredAbs(Scalar, Phi, Theta, Name=''):

        x, y, z = Sp2Cart(Scalar, Phi, Theta)

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


def PlotStructuredAmplitude(Scalar, Phi, Theta, Name=''):

    x, y, z = Sp2Cart(Phi*0+1, Phi, Theta)

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




def PlotField(Theta, Phi, Parallel, Perpendicular):

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10,6))
    axes[0,0].set_title('Parallel fields')
    axes[0,1].set_title('Perpendicular fields')
    axes[0,0].set_ylabel('Real Part')
    axes[1,0].set_ylabel('Imaginary Part')

    im0 = axes[0,0].pcolormesh(Theta,
                               Phi,
                               np.real(Parallel),
                               shading='auto')

    divider = make_axes_locatable(axes[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im0, cax=cax, orientation='vertical')

    im1 = axes[0,1].pcolormesh(Theta,
                               Phi,
                               np.real(Perpendicular),
                               shading='auto')

    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')


    im2 = axes[1,0].pcolormesh(Theta,
                               Phi,
                               np.imag(Parallel),
                               shading='auto')

    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im2, cax=cax, orientation='vertical')

    im3 = axes[1,1].pcolormesh(Theta,
                               Phi,
                               np.imag(Perpendicular),
                               shading='auto')

    divider = make_axes_locatable(axes[1,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im3, cax=cax, orientation='vertical')

    fig.tight_layout()
    plt.show(block=False)

# -
