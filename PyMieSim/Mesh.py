import numpy as np
from ai import cs
import math
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

from PyMieSim.Physics import Angle


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
                 MaxAngle:    float = np.pi/6,
                 Sampling:    int   = 1000,
                 PhiOffset:   float = 0,
                 GammaOffset: float = 0):

        self.PhiOffset = PhiOffset
        self.GammaOffset = GammaOffset
        self.MaxAngle = MaxAngle
        self.Sampling = Sampling
        Theta, Phi, dOmega = self.GenerateLedevedMesh()


        self.MakeProperties(Theta, Phi, dOmega)


    def MakeProperties(self, Theta, Phi, dOmega):

        self.Theta = Angle(Theta, unit='Radian')

        self.Phi = Angle(Phi, unit='Radian')

        self.Sampling = Theta.size

        self.SinMesh = np.abs( np.sin(Phi - np.pi/2) )

        self.dOmega = Angle(0, unit='Radian');
        self.dOmega.Radian = dOmega
        self.dOmega.Degree = dOmega * (180/np.pi)**2

        self.Omega = Angle(0, unit='Radian');
        self.Omega.Radian = np.abs( self.dOmega.Radian * self.Sampling)
        self.Omega.Degree = self.Omega.Radian * (180/np.pi)**2






    def Plot(self):
        x, y, z = cs.sp2cart(np.ones(self.Phi.Radian.shape),
                             self.Phi.Radian,
                             self.Theta.Radian)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X-direction [u.a.]')
        ax.set_ylabel('Y-direction [u.a.]')
        ax.set_zlabel('Z-direction [u.a.]')
        ax.scatter(x, y, z,s=10, c='k')
        ax.set_aspect('auto')
        ax.set_xlim([-1.3,1.3])
        ax.set_ylim([-1.3,1.3])
        ax.set_zlim([-1,1])

        ax.quiver(0,0,-2,0,0,1,length=0.5, color='k',arrow_length_ratio=0.5)
        ax.quiver(0,0,0,0,0,1.,length=1.5, color='k', arrow_length_ratio=0.1)
        ax.quiver(0,0,0,0,1,0,length=1.5, color='k', arrow_length_ratio=0.1)
        ax.quiver(0,0,0,1,0,0,length=1.5, color='k', arrow_length_ratio=0.1)

        phi, theta = np.mgrid[0.0:np.pi:20j, 0.0:2.0*np.pi:20j]

        ax.plot_surface(X          = np.sin(phi)*np.cos(theta),
                        Y          = np.sin(phi)*np.sin(theta),
                        Z          = np.cos(phi),
                        rstride    = 1,
                        cstride    = 1,
                        color      = 'b',
                        alpha      = 0.3,
                        linewidth  = 0.2,
                        shade      = True,
                        edgecolors = 'k'
                        )
        plt.show()


    def GenerateLedevedMesh(self):

        self.TrueSample, dOmega = self.ComputeTrueSample(self.Sampling)

        base = fibonacci_sphere(samples=self.Sampling, maxAngle=self.MaxAngle)

        base = self.AvoidPoles(base)

        r, phi, theta = cs.cart2sp(*base)

        base = cs.sp2cart(np.ones(phi.size), phi, theta)

        base = self.RotateOnPhi(-90, base)

        #base = self.RotateOnGamma(90 , base)

        r, phi, theta = cs.cart2sp(*base)

        self.base = Namespace(Phi=phi, Theta=theta)

        notbase = self.RotateOnPhi(self.PhiOffset, base)

        notbase = self.RotateOnGamma(self.GammaOffset, notbase)

        r, phi, theta = cs.cart2sp(*notbase)

        return theta, phi, dOmega



    def UpdateSphere(self, **kwargs):

        if 'MaxAngle' in kwargs: self.MaxAngle = kwargs['MaxAngle']
        if 'GammaOffset' in kwargs: self.GammaOffset = kwargs['GammaOffset']
        if 'PhiOffset' in kwargs: self.PhiOffset = kwargs['PhiOffset']
        if 'Sampling' in kwargs: self.Sampling = kwargs['Sampling']

        Theta, Phi, dOmega =  self.GenerateLedevedMesh()

        self.MakeProperties(Theta, Phi, dOmega)



    def AvoidPoles(self, base):
        Tinitial = cs.mx_rot_x(gamma = np.pi)

        return cs.mx_apply(Tinitial, *base)


    def RotateOnGamma(self, rotation, base):
        TPhi = cs.mx_rot_y(theta = rotation/180*np.pi)

        return cs.mx_apply(TPhi, *base)


    def RotateOnPhi(self, rotation, base):
        TGamma = cs.mx_rot_x(gamma = rotation/180*np.pi)

        return cs.mx_apply(TGamma, *base)


    def ComputeTrueSample(self, Sampling):

        SolidAngle = np.abs( 2*np.pi * (np.cos(self.MaxAngle) - np.cos(0)))

        ratio = (4*np.pi/SolidAngle)

        TrueSampling = int(Sampling*ratio)

        return TrueSampling, 4*np.pi / TrueSampling




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

        Phi, Theta = np.mgrid[-np.pi/2:np.pi/2:complex(Num), 0:2*np.pi:complex(Num)]

        self.Phi = Angle(Phi, unit='Radian')

        self.Theta = Angle(Theta, unit='Radian')






class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def fibonacci_sphere(samples=1, maxAngle=np.pi/2):

    X = []; Y = []; Z = []

    phi = math.pi * (3. - math.sqrt(5.))  ## golden angle in radians
    MaxY = np.cos(maxAngle)

    SolidAngle = np.abs( 2*np.pi * (np.cos(maxAngle) - np.cos(0)))

    ratio = 4*np.pi / SolidAngle

    TrueSampling = int( samples * ratio )

    for i in range(TrueSampling):
        y = 1 - (i / float(TrueSampling - 1)) * 2  ## y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  ## radius at y

        theta = phi * i  ## golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        X.append(x); Y.append(y); Z.append(z)

        if i >= samples - 1: break
        #if y <= MaxY: break

    return np.asarray(X), np.asarray(Y), np.asarray(Z)




# -
