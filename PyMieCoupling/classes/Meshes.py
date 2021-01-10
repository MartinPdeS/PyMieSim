import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from PyMieCoupling.functions.converts import rad2deg, deg2rad, NA2Angle
from PyMieCoupling.utils import Angle
from ai import cs
import math




class AngleMeshes(object):
    def __init__(self,
                 MaxAngle:    float = np.pi/6,
                 Samples:     int   = 1000,
                 PhiOffset          = 0,
                 GammaOffset        = 0):

        self.PhiOffset = PhiOffset
        self.GammaOffset = GammaOffset
        self.MaxAngle = MaxAngle
        self.Samples = Samples
        Theta, Phi, dOmega = self.GenerateLedevedMesh()


        self.MakeProperties(Theta, Phi, dOmega)


    def MakeProperties(self, Theta, Phi, dOmega):

        self.Theta = Angle(Theta, unit='Radian')

        self.Phi = Angle(Phi, unit='Radian')

        self.SinMesh = np.abs( np.sin(Phi - np.pi/2) )

        self.dOmega = Angle(0, unit='Radian')

        self.dOmega.Radian = dOmega

        self.dOmega.Degree = dOmega * (180/np.pi)**2


    def Plot(self):
        x, y, z = cs.sp2cart(np.ones(self.Phi.Radian.shape),
                             self.Phi.Radian,
                             self.Theta.Radian)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X-direction [u.a.]')
        ax.set_ylabel('Y-direction [u.a.]')
        ax.set_zlabel('Z-direction [u.a.]')
        ax.scatter(x, y, z,s=10,c='k')
        ax.set_aspect('auto')
        ax.set_xlim([-1.3,1.3])
        ax.set_ylim([-1.3,1.3])
        ax.set_zlim([-1,1])

        ax.quiver(-2,0,0,1,0,0,length=0.5, color='k',arrow_length_ratio=0.5)


        ax.quiver(0,0,0,0,0,1.,length=1.5, color='k', arrow_length_ratio=0.1)
        ax.quiver(0,0,0,0,1,0,length=1.5, color='k', arrow_length_ratio=0.1)
        ax.quiver(0,0,0,1,0,0,length=1.5, color='k', arrow_length_ratio=0.1)

        phi, theta = np.mgrid[0.0:np.pi:20j, 0.0:2.0*np.pi:20j]
        x = 1*np.sin(phi)*np.cos(theta)
        y = 1*np.sin(phi)*np.sin(theta)
        z = 1*np.cos(phi)
        ax.plot_surface(x,
                        y,
                        z,
                        rstride=1,
                        cstride=1,
                        color='b',
                        alpha=0.3,
                        linewidth=0.2,
                        shade=True,
                        edgecolors='k'
                        )


        plt.show()


    def GenerateLedevedMesh(self):

        #assert MaxAngle <= np.pi/2, print("Angle should be inferior to pi/2")

        self.TrueSample, dOmega = self.ComputeTrueSample(self.Samples)

        base = fibonacci_sphere(samples=self.Samples, maxAngle=self.MaxAngle)

        base = self.AvoidPoles(base)

        r, phi, theta = cs.cart2sp(*base)

        base = cs.sp2cart(np.ones(phi.size), phi, theta)

        base = self.RotateOnPhi(-90, base)

        base = self.RotateOnGamma(90 , base)

        r, phi, theta = cs.cart2sp(*base)

        self.base = Namespace(Phi=phi, Theta=theta) ###############################

        notbase = self.RotateOnPhi(self.PhiOffset, base)

        notbase = self.RotateOnGamma(self.GammaOffset, notbase)

        r, phi, theta = cs.cart2sp(*notbase)

        return theta, phi, dOmega


    def UpdateSphere(self, **kwargs):

        if 'MaxAngle' in kwargs: self.MaxAngle = kwargs['MaxAngle']
        if 'GammaOffset' in kwargs: self.GammaOffset = kwargs['GammaOffset']
        if 'PhiOffset' in kwargs: self.PhiOffset = kwargs['PhiOffset']
        if 'Samples' in kwargs: self.Samples = kwargs['Samples']

        Theta, Phi, dOmega =  self.GenerateLedevedMesh()

        self.MakeProperties(Theta, Phi, dOmega)


    def AvoidPoles(self, base):
        Tinitial = cs.mx_rot_x(gamma = np.pi)

        return cs.mx_apply(Tinitial, *base)


    def RotateOnGamma(self, rotation, base):
        TPhi = cs.mx_rot_z(phi = rotation/180*np.pi)

        return cs.mx_apply(TPhi, *base)


    def RotateOnPhi(self, rotation, base):
        TGamma = cs.mx_rot_x(gamma = rotation/180*np.pi)

        return cs.mx_apply(TGamma, *base)


    def ComputeTrueSample(self, Samples):

        SolidAngle = np.abs( 2*np.pi * (np.cos(self.MaxAngle) - np.cos(0)))

        ratio = (4*np.pi/SolidAngle)

        TrueSamples = int(Samples*ratio)

        return TrueSamples, 4*np.pi / TrueSamples


    def Cutoff(self, phi, theta):

        #indices = phi>=(np.pi/2-self.MaxAngle)

        #phi = phi[indices]; theta = theta[indices];

        #notbase = (phi*0+1, phi, theta)

        return cs.sp2cart(np.ones(phi.size), phi, theta)





class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)





def fibonacci_sphere(samples=1, maxAngle=np.pi/2):

    X = []; Y = []; Z = []

    phi = math.pi * (3. - math.sqrt(5.))  ## golden angle in radians
    MaxY = np.cos(maxAngle)

    solidAngle = np.pi * np.sin(maxAngle)**2

    ratio = 4*np.pi / solidAngle

    TrueSamples = int( samples * ratio )


    for i in range(TrueSamples):
        y = 1 - (i / float(TrueSamples - 1)) * 2  ## y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  ## radius at y

        theta = phi * i  ## golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        X.append(x); Y.append(y); Z.append(z)

        if y <= MaxY: break

    return np.asarray(X), np.asarray(Y), np.asarray(Z)



class DirectMeshes(object):

    def __init__(self,
                 XBound:   list = [-1,1],
                 YBound:   list = [-1,1],
                 XNpts:    int  = None,
                 YNpts:    int  = None,
                 Npts:     int  = None,
                 ):

        self.Npts = Npts
        self.XNpts = XNpts
        self.YNpts   = YNpts

        self.MakeProperties(XBound= XBound, YBound= YBound)


    def MakeProperties(self, XBound, YBound):

        XRange              = np.abs(XBound[0] - XBound[1])
        YRange                = np.abs(YBound[0] - YBound[1])

        XVector             = np.linspace(XBound[0], XBound[1], self.XNpts)
        YVector               = np.linspace(*YBound, self.YNpts)

        XMesh, YMesh      = np.meshgrid(XVector, YVector)

        XDelta = np.abs(XBound[0] - XBound[1])

        YDelta = np.abs(YBound[0] - YBound[1])

        self.dA   = XDelta * YDelta

        self.X = Namespace(Vector   =  XVector ,
                           Mesh     =  XMesh ,
                           Delta    =  XDelta ,
                            )

        self.Y = Namespace(Vector   =  YVector ,
                           Mesh     =  YMesh ,
                           Delta    =  YDelta
                              )




# -
