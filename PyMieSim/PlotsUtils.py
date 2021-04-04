import numpy as np
from numpy import pi, cos, sin, abs
from mayavi import mlab
from tvtk.tools import visual
from tvtk.api import tvtk
from tvtk.common import configure_input_data

from PyMieSim.utils import Sp2Cart

CMAP   = 'jet'
RED    = (1,0,0)
BLUE   = (0,0,1)
YELLOW = (1,1,0)
WHITE  = (1,1,1)

def ArrowAB(x1, y1, z1, x2, y2, z2, Scale=1):

    ar1             = visual.arrow(x = x1,
                                   y = y1,
                                   z = z1)
    ar1.length_cone = 0.2*Scale
    arrow_length    = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.actor.scale = [arrow_length*Scale, arrow_length*Scale, arrow_length*Scale]
    ar1.pos         = ar1.pos/arrow_length
    ar1.axis        = [x2-x1, y2-y1, z2-z1]


def ArrowAVec(Origin, Vec, Scale=0.5, Color=WHITE, ScaleTube=1.0):

    Vec = np.asarray(Vec)

    ar1  = visual.arrow(x     = Origin[0],
                        y     = Origin[1],
                        z     = Origin[2],
                        color = Color)

    ar1.length_cone = min( 1, 0.35/Scale*ScaleTube )

    ar1.actor.scale = np.array([1, 1, 1])*Scale

    ar1.radius_shaft = 0.02/Scale*ScaleTube

    ar1.radius_cone  = 2.5* ar1.radius_shaft

    ar1.pos         = np.asarray(Origin)/Scale

    ar1.axis        = Vec



def PlotCone(Figure,
             Radius      = 1.,
             Origin      = (0,0,0),
             Resolution  = 10,
             Height      = 1.,
             Direction   = (0,0,1)):

    cone = tvtk.ConeSource(center     = Origin,
                           radius     = Radius,
                           resolution = Resolution,
                           height     = Height,
                           direction  = Direction)

    configure_input_data( tvtk.PolyDataMapper(), cone.output )

    cone.update()

    p      = tvtk.Property(opacity=0.8, color=RED)

    actor  = tvtk.Actor(mapper=mapper, property=p)

    Figure.scene.add_actor(actor)


def ElectricFieldArrow(Figure, Origin, Pol, Scale=0.5):
    Vec = [sin(Pol),cos(Pol),0]

    ArrowAVec(Origin = Origin,
              Vec    = Vec,
              Scale  = Scale,
              Color  = RED)

    mlab.text3d(x          = Vec[0]+Origin[0],
                y          = Vec[1]+Origin[1],
                z          = Vec[2]+Origin[2],
                text       = 'E',
                line_width = 0.1,
                figure     = Figure,
                scale      = 0.25,
                color      = RED)

def WavenumberArrow(Figure, Origin, Scale=0.5):
    Vec = [0,0,1]

    ArrowAVec(Origin = Origin,
              Vec    = Vec,
              Scale  = Scale,
              Color  = YELLOW)

    mlab.text3d(x          = Vec[0]+Origin[0],
                y          = Vec[1]+Origin[1],
                z          = Vec[2]+Origin[2],
                text       = 'k',
                line_width = 0.1,
                figure     = Figure,
                scale      = 0.25,
                color      = YELLOW)


def MagneticFieldArrow(Figure, Origin, Pol, Scale=0.5):
    Vec = [cos(Pol),sin(Pol),0]

    ArrowAVec(Origin = Origin,
              Vec    = Vec,
              Scale  = Scale,
              Color  = BLUE)

    mlab.text3d(x          = Vec[0]+Origin[0],
                y          = Vec[1]+Origin[1],
                z          = Vec[2]+Origin[2],
                text       = 'B',
                line_width = 0.1,
                figure     = Figure,
                scale      = 0.25,
                color      = BLUE)


def AddSource(Figure, Origin, Polarization, Scale=0.5):
    ElectricFieldArrow(Figure, Origin, Polarization, Scale=Scale)
    MagneticFieldArrow(Figure, Origin, Polarization, Scale=Scale)
    WavenumberArrow(Figure, Origin, Scale=Scale)


def AddUnitSphere(Num, Radius, Origin, Figure):

    Coord = list( UnitSphere(Num=50, Radius=Radius) )

    Coord[0] += Origin[0]
    Coord[1] += Origin[1]
    Coord[2] += Origin[2]

    mlab.mesh(*Coord, colormap='gray', opacity=0.5, figure=Figure)


def AddUnitAxes(Figure,
                Scale      = 1.,
                Origin     = (0,0,0),
                Label      = True,
                ScaleTube  = 1.0):

    mlab.points3d(0, 0, 0, opacity=0.01, scale_factor=0.01, figure=Figure)

    ArrowAVec(Origin, (0, 0, 1), Scale=Scale, ScaleTube=ScaleTube)
    if Label: mlab.text3d(x          = Origin[0],
                          y          = Origin[1],
                          z          = Origin[2]+Scale,
                          text       = 'Z',
                          line_width = 0.1,
                          figure     = Figure,
                          scale      = 0.25)

    ArrowAVec(Origin, (0, 1, 0), Scale=Scale, ScaleTube=ScaleTube)
    if Label: mlab.text3d(x          = Origin[0],
                          y          = Origin[1]+Scale,
                          z          = Origin[2],
                          text       = 'Y',
                          line_width = 0.1,
                          figure     = Figure,
                          scale      = 0.25)

    ArrowAVec(Origin, (1, 0, 0), Scale=Scale, ScaleTube=ScaleTube)
    if Label: mlab.text3d(x          = Origin[0]+Scale,
                          y          = Origin[1],
                          z          = Origin[2],
                          text       = 'X',
                          line_width = 0.01,
                          figure     = Figure,
                          scale      = 0.25)


def UnitSphere(Num, Radius=1.):
    phi, theta = np.mgrid[-pi/2:pi/2:complex(Num), -pi:pi:complex(Num)]

    return Sp2Cart(phi*0+Radius, phi, theta)


def implicit_plot(expr, ext_grid, fig_handle=None, Nx=101, Ny=101, Nz=101,
                 col_isurf=(50/255, 199/255, 152/255), col_osurf=(240/255,36/255,87/255),
                 opa_val=0.8, opaque=True, ori_axis=True, **kwargs):
    """Function to plot algebraic surfaces described by implicit equations in Mayavi

    Implicit functions are functions of the form

        `F(x,y,z) = c`

    where `c` is an arbitrary constant.

    Parameters
    ----------
    expr : string
        The expression `F(x,y,z) - c`; e.g. to plot a unit sphere, the `expr` will be `x**2 + y**2 + z**2 - 1`
    ext_grid : 6-tuple
        Tuple denoting the range of `x`, `y` and `z` for grid; it has the form - (xmin, xmax, ymin, ymax, zmin, zmax)
    fig_handle : figure handle (optional)
        If a mayavi figure object is passed, then the surface shall be added to the scene in the given figure. Then, it is the responsibility of the calling function to call mlab.show().
    Nx, Ny, Nz : Integers (optional, preferably odd integers)
        Number of points along each axis. It is recommended to use odd numbers to ensure the calculation of the function at the origin.
    col_isurf : 3-tuple (optional)
        color of inner surface, when double-layered surface is used. This is also the specified color for single-layered surface.
    col_osurf : 3-tuple (optional)
        color of outer surface
    opa_val : float (optional)
        Opacity value (alpha) to use for surface
    opaque : boolean (optional)
        Flag to specify whether the surface should be opaque or not
    ori_axis : boolean
        Flag to specify whether a central axis to draw or not

    """
    if fig_handle==None:  # create a new figure
        fig = mlab.figure(1,bgcolor=(0.97, 0.97, 0.97), fgcolor=(0, 0, 0), size=(800, 800))
    else:
        fig = fig_handle
    xl, xr, yl, yr, zl, zr = ext_grid
    x, y, z = np.mgrid[xl:xr:eval('{}j'.format(Nx)),
                       yl:yr:eval('{}j'.format(Ny)),
                       zl:zr:eval('{}j'.format(Nz))]
    scalars = eval(expr)
    src = mlab.pipeline.scalar_field(x, y, z, scalars)
    if opaque:
        delta = 1.e-5
        opa_val=1.0
    else:
        delta = 0.0
        #col_isurf = col_osurf
    # In order to render different colors to the two sides of the algebraic surface,
    # the function plots two contour3d surfaces at a &quot;distance&quot; of delta from the value
    # of the solution.
    # the second surface (contour3d) is only drawn if the algebraic surface is specified
    # to be opaque.
    cont1 = mlab.pipeline.iso_surface(src, color=col_isurf, contours=[0-delta],
                                      transparent=False, opacity=opa_val)
    cont1.compute_normals = False # for some reasons, setting this to true actually cause
                                  # more unevenness on the surface, instead of more smooth
    if opaque: # the outer surface is specular, the inner surface is not
        cont2 = mlab.pipeline.iso_surface(src, color=col_osurf, contours=[0+delta],
                                          transparent=False, opacity=opa_val)
        cont2.compute_normals = False
        cont1.actor.property.backface_culling = True
        cont2.actor.property.frontface_culling = True
        cont2.actor.property.specular = 0.2 #0.4 #0.8
        cont2.actor.property.specular_power = 55.0 #15.0
    else:  # make the surface (the only surface) specular
        cont1.actor.property.specular = 0.2 #0.4 #0.8
        cont1.actor.property.specular_power = 55.0 #15.0

    # Scene lights (4 lights are used)
    engine = mlab.get_engine()
    scene = engine.current_scene
    cam_light_azimuth = [78, -57, 0, 0]
    cam_light_elevation = [8, 8, 40, -60]
    cam_light_intensity = [0.72, 0.48, 0.60, 0.20]
    """
    for i in range(4):
        camlight = scene.scene.light_manager.lights[i]
        camlight.activate = True
        camlight.azimuth = cam_light_azimuth[i]
        camlight.elevation = cam_light_elevation[i]
        camlight.intensity = cam_light_intensity[i].
    """
    # axis through the origin
    if ori_axis:
        len_caxis = int(1.05*np.max(np.abs(np.array(ext_grid))))
        caxis = mlab.points3d(0.0, 0.0, 0.0, len_caxis, mode='axes',color=(0.15,0.15,0.15),
                              line_width=1.0, scale_factor=1.,opacity=1.0)
        caxis.actor.property.lighting = False
    # if no figure is passed, the function will create a figure.
    if fig_handle==None:
        # Setting camera
        cam = fig.scene.camera
        cam.elevation(-20)
        cam.zoom(1.0) # zoom should always be in the end.
        mlab.show()


# -
