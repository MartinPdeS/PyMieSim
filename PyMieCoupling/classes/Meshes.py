import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from PyMieCoupling.functions.converts import rad2deg, deg2rad, NA2Angle
from PyMieCoupling.utils import Angle
from ai import cs
import math





def fibonacci_sphere(samples=1):

    X = []; Y = []; Z = []

    phi = math.pi * (3. - math.sqrt(5.))  ## golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  ## y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  ## radius at y

        theta = phi * i  ## golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        X.append(x)
        Y.append(y)
        Z.append(z)

    return np.asarray(X), np.asarray(Y), np.asarray(Z)




class AngleMeshes(object):
    def __init__(self,
                 MaxAngle:    float = np.pi/6,
                 Samples:     int   = 1000,
                 PhiOffset          = 0,
                 GammaOffset        = 0):

        theta, phi, dOmega = self.GenerateLedevedMesh(MaxAngle    = MaxAngle,
                                                      Samples     = Samples,
                                                      PhiOffset   = PhiOffset,
                                                      GammaOffset = GammaOffset)


        self.Theta = Angle(theta, unit='Radian')

        self.Phi = Angle(phi, unit='Radian')

        self.SinMesh = np.sin(phi)

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

        ax.quiver(0,0,-1.5,0,0,1,length=0.5, color='k')

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


    def GenerateLedevedMesh(self, MaxAngle, Samples, PhiOffset, GammaOffset):

        #assert MaxAngle <= np.pi/2, print("Angle should be inferior to pi/2")

        TrueSample, dOmega = self.ComputeTrueSample(Samples, MaxAngle)

        base = fibonacci_sphere(samples=TrueSample)

        Tinitial = cs.mx_rot_x(gamma = 90/180*np.pi)

        Tfinal = cs.mx_rot_x(gamma = -90/180*np.pi)

        base = cs.mx_apply(Tinitial, *base)

        base = self.Cutoff(*base, MaxAngle)

        notbase = self.RotateOnPhi(PhiOffset, base)

        notbase = self.RotateOnTheta(GammaOffset, notbase)

        #notbase = cs.mx_apply(Tfinal, *notbase)

        r, phi, theta = cs.cart2sp(*notbase)

        return theta, phi, dOmega


    def RotateOnTheta(self, rotation, base):
        TPhi = cs.mx_rot_y(theta = rotation/180*np.pi)

        return cs.mx_apply(TPhi, *base)


    def RotateOnPhi(self, rotation, base):
        TTheta = cs.mx_rot_x(gamma = rotation/180*np.pi)

        return cs.mx_apply(TTheta, *base)


    def ComputeTrueSample(self, Samples, MaxAngle):

        SolidAngle = np.abs( 2*np.pi * (np.cos(MaxAngle) - np.cos(0)))

        ratio = (4*np.pi/SolidAngle)

        TrueSamples = int(Samples*ratio)

        return TrueSamples, 4*np.pi / TrueSamples


    def Cutoff(self, x, y, z, MaxAngle):
        T0 = cs.mx_rot_x(gamma = 90/180*np.pi)

        x,y,z = cs.mx_apply(T0, x, y, z)

        r, phi, theta = cs.cart2sp(x, y, z)

        indices = phi>=(np.pi/2-MaxAngle)

        phi = phi[indices]; theta = theta[indices]; r = r[indices]

        self.base = (r,phi,theta)

        return cs.sp2cart(r, phi, theta)








class _AngleMeshes(object):

    def __init__(self,
                 ThetaBound:  list  = [-180,180],
                 PhiBound:    list  = [-180,180],
                 ThetaNpts:   int   = None,
                 PhiNpts:     int   = None,
                 GammaOffset: float = 0,
                 PhiOffset:   float = 0
                 ):

        self.ThetaNpts = ThetaNpts
        self.PhiNpts   = PhiNpts
        self.MakeProperties(ThetaBound  = ThetaBound,
                            PhiBound    = PhiBound,
                            GammaOffset = GammaOffset,
                            PhiOffset   = PhiOffset)


    def MakeProperties(self, ThetaBound, PhiBound, GammaOffset, PhiOffset):

        ThetaMesh, PhiMesh = np.mgrid[ThetaBound[0] : ThetaBound[1] : complex(0,self.ThetaNpts),
                                      PhiBound[0]   : PhiBound[1]   : complex(self.PhiNpts), ]

        ThetaDelta = np.abs(ThetaBound[0] - ThetaBound[1])

        PhiDelta = np.abs(PhiBound[0] - PhiBound[1])

        self.dOmega = Angle( 0 )

        self.dOmega.Degree = PhiDelta * ThetaDelta

        self.dOmega.Radian = deg2rad(PhiDelta) * deg2rad(ThetaDelta)

        self.Theta = Namespace( Boundary = Angle( [ThetaBound[0], ThetaBound[1]] ),
                                Offset   = Angle( GammaOffset ),
                                Mesh     = Angle( ThetaMesh ),
                                Delta    = Angle( ThetaDelta ),
                                )

        self.Phi = Namespace( Boundary = Angle( [PhiBound[0], PhiBound[1]] ),
                              Offset   = Angle( PhiOffset ),
                              Mesh     = Angle( PhiMesh ),
                              Delta    = Angle( PhiDelta )
                              )

        self.SinMesh = np.abs(np.sin(PhiMesh))

        if PhiOffset != 0: self.MakePhiOffset(PhiOffset)


    def MakePhiOffset(self, rotation):

        self.PhiOffset = Angle(rotation)

        r = np.ones(self.Phi.Mesh.Radian.flatten().shape)


        x, y, z = cs.sp2cart(r,
                              self.Phi.Mesh.Radian.flatten()-np.pi/2,
                              self.Theta.Mesh.Radian.flatten(),
                              )

        Tz = cs.mx_rot_x(gamma = rotation/180*np.pi)

        xp,yp,zp = cs.mx_apply(Tz, x, y, z)

        rp, phip, thetap = cs.cart2sp(xp, yp, zp)

        self.Theta.Mesh   = Angle( thetap.reshape(self.Theta.Mesh.Radian.shape), unit='Radian' )

        self.Phi.Mesh = Angle( (phip.reshape(self.Phi.Mesh.Radian.shape)+np.pi/2) , unit='Radian' )





    def Plot(self):

        R = np.ones(self.Theta.Mesh.Radian.flatten().shape)
        xp, yp, zp = cs.sp2cart(R,
                                self.Phi.Mesh.Radian.flatten()+np.pi/2,
                                self.Theta.Mesh.Radian.flatten(),
                                )


        r = 1
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:50j, 0.0:2.0*pi:50j]
        x = r*sin(phi)*cos(theta)
        y = r*sin(phi)*sin(theta)
        z = r*cos(phi)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(x,
                        y,
                        z,
                        rstride=1,
                        cstride=1,
                        color='k',
                        alpha=0.3,
                        linewidth=0.00,
                        edgecolors='k',
                        antialiased = False)

        ax.scatter(xp,yp,zp,color="k",s=2)
        ax.set_xlabel('X-Direction')
        ax.set_ylabel('Y-Direction')
        ax.set_zlabel('Z-Direction')
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
        ax.set_xlim([-1,1]); ax.set_ylim([-1,1]); ax.set_zlim([-1,1])
        plt.tight_layout()
        plt.show()



class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)



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
