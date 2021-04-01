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


def PlotUnitSphere(Num, Radius, Origin, Figure):
    Coord = UnitSphere(Num=50, Radius=1.)

    Coord = list(Coord)

    Coord[0] += Origin[0]
    Coord[1] += Origin[1]
    Coord[2] += Origin[2]

    mlab.mesh(*Coord, colormap='gray', opacity=0.5, figure=Figure)


def ArrowAVec(Origin, Vec, Scale=0.5, Color=(1,1,1)):

    ar1             = visual.arrow(x     = Origin[0],
                                   y     = Origin[1],
                                   z     = Origin[2],
                                   color = Color)

    ar1.length_cone = 0.35*Scale
    ar1.actor.scale = np.array([1, 1, 1])* Scale
    ar1.pos         = np.asarray(Origin)/Scale
    ar1.axis        = Vec


def PlotUnitAxes(Figure, Scale=1., Origin=(0,0,0)):

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


def Unstructured(**kwargs):
    if kwargs['Mode'] == 'Absolute':
        kwargs.pop('Mode')
        return UnstructuredAbs(**kwargs)

    if kwargs['Mode'] == 'Amplitude':
        kwargs.pop('Mode')
        return UnstructuredAmplitude(**kwargs)


def Structured(**kwargs):
    if kwargs['Mode'] == 'Absolute':
        kwargs.pop('Mode')
        return StructuredAbs(**kwargs)

    if kwargs['Mode'] == 'Amplitude':
        kwargs.pop('Mode')
        return StructuredAmplitude(**kwargs)

@mlab.show
def UnstructuredAbs(Mesh, Scalar=None, Name=''):

    if Scalar is None: Scalar = Mesh.Phi.Radian*0+1

    x, y, z = Sp2Cart(Scalar.flatten(),
                      Mesh.Phi.Radian.flatten(),
                      Mesh.Theta.Radian.flatten())

    fig = mlab.figure(figure = Name, size=(600,300))
    visual.set_viewer(fig)

    PlotUnitAxes(Figure=fig, Scale=1., Origin=(0,0,-1.5))

    PlotUnitSphere(Num=50, Radius=1, Origin=(0,0,0), Figure=fig)

    im0 = mlab.points3d(x,
                        y,
                        z,
                        abs(Scalar),
                        mode       = 'sphere',
                        scale_mode = 'none',
                        colormap   = CMAP)

    mlab.colorbar(object      = im0,
                  label_fmt   = "%.0e",
                  nb_labels   = 5,
                  title       = 'Real part',
                  orientation = 'horizontal' )


@mlab.show
def UnstructuredAmplitude(Mesh, Scalar=None, Name=''):

    if Scalar is None: Scalar = Mesh.Phi.Radian.flatten()*0+1

    x, y, z = Sp2Cart(Scalar,
                      Mesh.Phi.Radian.flatten(),
                      Mesh.Theta.Radian.flatten())

    fig = mlab.figure(figure=Name,size=(600,300))
    visual.set_viewer(fig)

    O0 = (+3, 0, 0)
    O1 = (-3, 0, 0)

    PlotUnitAxes(Figure=fig, Scale=1., Origin=(+3, 0, -2))
    PlotUnitAxes(Figure=fig, Scale=1., Origin=(-3, 0, -2))

    mlab.text(O0[0], O0[1],  u'Imaginary', z=3., width=0.15)
    mlab.text(O1[0], O1[1], u'Real', z=3,  width=0.07)

    PlotUnitSphere(Num    = 50,
                   Radius = 1.,
                   Origin = O0,
                   Figure = fig)

    PlotUnitSphere(Num    = 50,
                   Radius = 1.,
                   Origin = O1,
                   Figure = fig)


    im0 = mlab.points3d(Mesh.CartCoord[0]+O0[0],
                        Mesh.CartCoord[1]+O0[1],
                        Mesh.CartCoord[2]+O0[2],
                        Scalar.real,
                        mode       = 'sphere',
                        scale_mode = 'none',
                        colormap   = CMAP)

    im1 = mlab.points3d(Mesh.CartCoord[0]+O1[0],
                        Mesh.CartCoord[1]+O1[1],
                        Mesh.CartCoord[2]+O1[2],
                        Scalar.imag,
                        mode        = 'sphere',
                        scale_mode  = 'none',
                        colormap    = CMAP)

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

@mlab.show
def StructuredAbs(Scalar, Phi, Theta, Name='', Polarization=None):

    fig = mlab.figure(figure=Name,size=(600,400))
    visual.set_viewer(fig)

    O = (0,0,0)

    if Polarization is not None: PolarizationArrow(O, Polarization, Scale=0.4)

    x, y, z = Sp2Cart(Scalar/np.max(Scalar), Phi, Theta)

    PlotUnitAxes(Figure=fig, Scale=0.3, Origin=O)

    im = mlab.mesh(x, y, z-np.min(z)+0.5, colormap='viridis', figure=fig)

    mlab.colorbar(object      = im,
                  label_fmt   = "%.0e",
                  nb_labels   = 5,
                  title       = 'Normalized scale',
                  orientation = 'horizontal' )


@mlab.show
def StructuredAmplitude(Scalar, Phi, Theta, Name='', Polarization=None):

    x, y, z = Sp2Cart(Phi*0+1, Phi, Theta)

    Xoffset = 3
    O0      = (+3,0,0)
    O1      = (-3,0,0)

    fig = mlab.figure(figure=Name,size=(600,300))
    visual.set_viewer(fig)

    if Polarization is not None:
        PolarizationArrow(O0, Polarization, Scale=1.1)
        PolarizationArrow(O1, Polarization, Scale=1.1)

    PlotUnitAxes(Figure=fig, Origin=O0, Scale=1)
    PlotUnitAxes(Figure=fig, Origin=O1, Scale=1)

    mlab.text(O0[0], O0[1], u'Imaginary', z=5., width=0.08)
    mlab.text(O1[0], O1[1], u'Real',      z=5., width=0.05)

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
