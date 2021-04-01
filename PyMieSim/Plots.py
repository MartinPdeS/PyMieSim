import numpy as np
from numpy import pi, cos, sin, abs
from mayavi import mlab
from tvtk.tools import visual

from PyMieSim.Physics import Angle
from PyMieSim.utils import Sp2Cart, Cart2Sp

CMAP = 'jet'

def ArrowAB(x1, y1, z1, x2, y2, z2, Scale=1):

    ar1             = visual.arrow(x = x1,
                                   y = y1,
                                   z = z1)
    ar1.length_cone = 0.2*Scale
    arrow_length    = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.actor.scale = [arrow_length*Scale, arrow_length*Scale, arrow_length*Scale]
    ar1.pos         = ar1.pos/arrow_length
    ar1.axis        = [x2-x1, y2-y1, z2-z1]



def PolarizationArrow(Origin, Pol, Scale=0.5):

    return ArrowAVec(Origin = Origin,
                     Vec    = [sin(Pol),cos(Pol),0],
                     Scale  = Scale,
                     Color  = (1,0.2,0.2))


def ArrowAVec(Origin, Vec, Scale=0.5, Color=(1,1,1)):

    ar1             = visual.arrow(x     = Origin[0],
                                   y     = Origin[1],
                                   z     = Origin[2],
                                   color = Color)
    ar1.length_cone = 0.35*Scale
    ar1.actor.scale = np.array([1, 1, 1])* Scale
    ar1.pos         = np.asarray(Origin)/Scale
    ar1.axis        = Vec


def UnitAxes(Figure, Scale=1., Origin=(0,0,0)):

    mlab.plot3d(0,0,0, line_width=1e-12, figure=Figure)

    ArrowAVec(Origin, (0, 0, 1), Scale=Scale )
    ArrowAVec(Origin, (0, 1, 0), Scale=Scale )
    ArrowAVec(Origin, (1, 0, 0), Scale=Scale )

    mlab.text(Origin[0], Origin[1], 'Z', z=Origin[2]+Scale, width=0.01, figure=Figure)
    mlab.text(Origin[0], Origin[1]+Scale, 'Y', z=Origin[2], width=0.01, figure=Figure)
    mlab.text(Origin[0]+Scale, Origin[1], 'X', z=Origin[2], width=0.01, figure=Figure)


def UnitSphere(Num, Radius=1.):

    phi, theta = np.mgrid[-pi/2:pi/2:complex(Num), -pi:pi:complex(Num)]
    return Sp2Cart(phi*0+Radius, phi, theta)


def PlotUnstructuredAbs(Scalar, Mesh, Name=''):

    x, y, z = Sp2Cart(Mesh.Phi.Radian.flatten()*0+1,
                      Mesh.Phi.Radian.flatten(),
                      Mesh.Theta.Radian.flatten())

    _PlotUnstructuredAbs(Scalar  = Mesh.Phi.Radian*0+1,
                         Phi     = Mesh.Phi.Radian,
                         Theta   = Mesh.Theta.Radian,
                         Name    = Name)


def _PlotUnstructuredAbs(Scalar, Phi, Theta, Name=''):

    Phi     = Phi.flatten()
    Theta   = Theta.flatten()
    Scalar  = Scalar.flatten()

    x, y, z = Sp2Cart(Phi*0+1,
                      Phi,
                      Theta)

    offset   = 3
    X        = np.linspace(-1,1,10)*1.3
    Y        = np.linspace(-1,1,10)*1.3
    Z        = np.linspace(-1,1,10)*2

    zeros    = np.zeros_like(X)

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

    im0 = mlab.points3d(x + offset, y, z, np.abs(Scalar), mode='sphere', scale_mode='none', colormap=CMAP)

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

    UnitAxes(Figure=fig, Scale=1., )

    UnitAxes(Figure=fig, Scale=1., Origin=(0, 0,-3))

    mlab.text(offset, 0, u'Imaginary', z=-2., width=0.15)
    mlab.text(0, 0, u'Real', z=-2, width=0.07)

    mlab.text(offset/2, 0, Name, z=2, width=0.2)

    xp, yp, zp = UnitSphere(Num=50, Radius=1.)

    mlab.mesh(xp, yp, zp, colormap='gray', opacity=0.5)

    mlab.mesh(xp+offset, yp, zp, colormap='gray', opacity=0.5)

    im0 = mlab.points3d(*Mesh.CartCoord, Scalar.real, mode='sphere', scale_mode='none', colormap=CMAP)

    im1 = mlab.points3d(Mesh.CartCoord[0]+offset, Mesh.CartCoord[1], Mesh.CartCoord[2], Scalar.imag, mode='sphere', scale_mode='none', colormap=CMAP)

    mlab.colorbar(object = im0, label_fmt="%.0e", nb_labels=5, title='Real part', orientation='horizontal' )

    mlab.colorbar(object = im1, label_fmt="%.0e", nb_labels=5, title='Imaginary part', orientation='vertical' )


def StructuredAbs(Scalar, Phi, Theta, Name='', Polarization=None):

    fig = mlab.figure(figure=Name,size=(600,600))
    visual.set_viewer(fig)

    O = (0,0,0)

    if Polarization: PolarizationArrow(O, Polarization, Scale=0.6)

    x, y, z = Sp2Cart(Scalar/np.max(Scalar), Phi, Theta)

    UnitAxes(Figure=fig, Scale=0.4, Origin=O)

    mlab.mesh(x, y, z+2, colormap='viridis', figure=fig)



def StructuredAmplitude(Scalar, Phi, Theta, Name='', Polarization=None):

    x, y, z = Sp2Cart(Phi*0+1, Phi, Theta)

    Xoffset = 3
    O0      = (Xoffset,0,0)
    O1      = (-Xoffset,0,0)

    fig = mlab.figure(figure=Name,size=(600,300))

    visual.set_viewer(fig)

    if Polarization:
        PolarizationArrow(O0, Polarization, Scale=1.1)
        PolarizationArrow(O1, Polarization, Scale=1.1)

    UnitAxes(Figure=fig, Origin=O0, Scale=1)

    UnitAxes(Figure=fig, Origin=O1, Scale=1)

    mlab.text(-Xoffset+1, 0, u'Imaginary', z=5., width=0.08)
    mlab.text(+Xoffset+1, 0, u'Real',      z=5., width=0.05)

    im0 = mlab.points3d(x-Xoffset,
                        y,
                        z+2,
                        Scalar.real,
                        mode       = 'sphere',
                        scale_mode = 'none',
                        colormap   = CMAP)

    im1 = mlab.points3d(x+Xoffset,
                        y,
                        z+2,
                        Scalar.imag,
                        mode       = 'sphere',
                        scale_mode = 'none',
                        colormap   = CMAP)

    mlab.colorbar(object      = im0,
                  label_fmt   = "%.0e",
                  nb_labels   = 5,
                  title       = 'Real part',
                  orientation = 'horizontal' )

    mlab.colorbar(object      = im1,
                  label_fmt   = "%.0e",
                  nb_labels   = 5,
                  title       = 'Imaginary part',
                  orientation = 'vertical' )



# -
