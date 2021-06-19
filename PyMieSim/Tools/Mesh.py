import numpy as np
from mayavi import mlab

from PyMieSim.Physics         import Angle
from PyMieSim.Tools.Plots     import Unstructured
from PyMieSim.Tools.Fibonacci import Mesh as FMesh
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


        self.Structured  = False
        self.Sampling    = Sampling
        self.MaxAngle    = MaxAngle
        self.PhiOffset   = PhiOffset
        self.GammaOffset = GammaOffset
        self.GenerateLedevedMesh()


    def MakeProperties(self):

        self.CartCoord = np.asarray([self.bind.x, self.bind.y, self.bind.z])

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

        Unstructured(Mesh=self, Name=Name, Mode='Amplitude')

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

        self.Structure     = True

        self.Num           = Num
        self.shape         = [Num, Num]

        Phi, Theta         = np.mgrid[-pi/2:pi/2:complex(Num), 0:2*pi:complex(Num)]

        self.dOmega         = Angle(0, unit='Radian');
        self.dOmega.Radian  = np.abs( Phi[0,0] - Phi[1,0] ) * np.abs( Theta[0,0] - Theta[0,1] )
        self.dOmega.Degree  = self.dOmega.Radian * (180/pi)**2


        self.Omega         = Angle(0, unit='Radian');
        self.Omega.Radian  = 4 * pi
        self.Omega.Degree  = self.dOmega.Radian * (180/pi)**2

        self.Phi           = Angle(Phi.flatten(), unit='Radian')
        self.Theta         = Angle(Theta.flatten(), unit='Radian')

        self.SinMesh       = np.sin(self.Phi.Radian+pi/2)




class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)



# -
