import numpy as np
from mayavi import mlab

from PyMieSim.Plots import UnitSphere, UnitAxes
from PyMieSim.Physics import Angle
from PyMieSim.Fibonacci import Mesh as FMesh
pi = np.pi



class FibonacciMesh(object):
    """
    Class wich represent an angular mesh. The distribution of points inside
    the mesh is similar to a Fibonacci sphere where each point cover an
    equivalent solid angle.

    Parameters
    ----------

    MaxAngle : float
        Angle in radian defined by the numerical aperture of the imaging system.
    Sampling : int
        Number of point distrubuted inside the Solid angle defined by the
        numerical aperture.
    PhiOffset : float
        Angle offset in the parallel direction of the polarization of
        incindent light.
    GammaOffset : float
        Angle offset in the perpendicular direction of the polarization of
        incindent light.

    """

    def __init__(self,
                 MaxAngle:    float = 1.5,
                 Sampling:    int   = 1000,
                 PhiOffset:   float = 0.,
                 GammaOffset: float = 0.):


        self.Sampling    = Sampling
        self.MaxAngle    = MaxAngle
        self.PhiOffset   = PhiOffset
        self.GammaOffset = GammaOffset
        self.GenerateLedevedMesh()


    def MakeProperties(self):

        self.CartCoord = (self.bind.x, self.bind.y, self.bind.z)

        self.SphCoord  = (self.bind.r, self.bind.phi, self.bind.theta)

        self.Theta         = self.bind.r
        self.Theta         = Angle(self.bind.theta, unit='Radian')
        self.Phi           = Angle(self.SphCoord[1], unit='Radian')

        self.SinMesh       = np.abs( np.sin( self.SphCoord[1] - np.pi/2 ) )

        self.dOmega        = Angle(0, unit='Radian');
        self.dOmega.Radian = self.bind.dOmega
        self.dOmega.Degree = self.bind.dOmega * (180/pi)**2

        self.Omega         = Angle(0, unit='Radian');
        self.Omega.Radian  = self.bind.Omega
        self.Omega.Degree  = self.bind.Omega * (180/pi)**2


    def Plot(self):

        Name = 'Angular mesh'

        fig = mlab.figure(size=(600,300))

        UnitAxes(fig)

        mlab.text(0, 0, Name, z=-2, width=0.2)

        coord = UnitSphere(Num=50, Radius=1.)

        mlab.mesh(*coord, colormap='gray', opacity=0.5)

        im0 = mlab.points3d(*self.CartCoord,
                             mode         = 'sphere',
                             scale_mode   = 'none',
                             scale_factor = 0.03)

        mlab.show()


    def GenerateLedevedMesh(self):

        self.bind = FMesh(self.Sampling,
                          self.MaxAngle,
                          np.deg2rad(self.PhiOffset),
                          np.deg2rad(self.GammaOffset) )


        self.base = (self.bind.PhiBase, self.bind.ThetaBase)

        self.MakeProperties()


    def UpdateSphere(self, **kwargs):

        if 'MaxAngle'    in kwargs: self.MaxAngle    = kwargs['MaxAngle']
        if 'GammaOffset' in kwargs: self.GammaOffset = kwargs['GammaOffset']
        if 'PhiOffset'   in kwargs: self.PhiOffset   = kwargs['PhiOffset']
        if 'Sampling'    in kwargs: self.Sampling    = kwargs['Sampling']

        self.GenerateLedevedMesh()


class StructuredFullMesh(object):
    """Class wich represent an angular mesh.

    Parameters
    ----------
    MaxAngle : float
        Angle in radian defined by the numerical aperture of the imaging system.
    Num : int
        Number of point for the mesh shape [Num, Num].

    """
    def __init__(self, Num: int = 100):
        self.Num = Num

        self.shape = [Num, Num]

        Phi, Theta = np.mgrid[-pi/2:pi/2:complex(Num), 0:2*pi:complex(Num)]

        self.Phi = Angle(Phi, unit='Radian')

        self.Theta = Angle(Theta, unit='Radian')




class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)



# -
