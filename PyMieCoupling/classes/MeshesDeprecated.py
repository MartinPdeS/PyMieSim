import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.functions.converts import rad2deg, deg2rad
from PyMieCoupling.utils import Angle
from ai import cs

class AngleMeshes(object):

    def __init__(self,
                 ThetaBound:  list = [-180,180],
                 PhiBound:    list = [-180,180],
                 ThetaNpts:   int  = None,
                 PhiNpts:     int  = None,
                 ThetaOffset: float = 0,
                 PhiOffset:   float = 0
                 ):

        self.ThetaNpts = ThetaNpts
        self.PhiNpts   = PhiNpts
        self.MakeProperties(ThetaBound  = ThetaBound,
                            PhiBound    = PhiBound,
                            ThetaOffset = ThetaOffset,
                            PhiOffset   = PhiOffset)


    def MakeProperties(self,
                       ThetaBound,
                       PhiBound,
                       ThetaOffset,
                       PhiOffset):


        if self.PhiNpts:
            ThetaVector             = np.linspace(ThetaBound[0], ThetaBound[1], self.ThetaNpts)
            PhiVector               = np.linspace(*PhiBound, self.PhiNpts)

        else:
            ThetaVector             = np.linspace(ThetaBound[0], ThetaBound[1], self.Npts)
            PhiVector               = np.linspace(*PhiBound, self.Npts)
        ThetaMesh, PhiMesh      = np.meshgrid(ThetaVector, PhiVector)

        ThetaDelta = np.abs(ThetaBound[0] - ThetaBound[1])

        PhiDelta = np.abs(PhiBound[0] - PhiBound[1])

        self.dOmega   = Angle( 0 )

        self.dOmega.Degree = PhiDelta * ThetaDelta

        self.dOmega.Radian = deg2rad(PhiDelta) * deg2rad(ThetaDelta)

        self.Theta = Namespace( Offset   = Angle( ThetaOffset ),
                                Vector   = Angle( ThetaVector ),
                                Mesh     = Angle( ThetaMesh ),
                                Delta    = Angle( ThetaDelta ),
                                )

        self.Phi = Namespace( Offset   = Angle( PhiOffset ),
                              Vector   = Angle( PhiVector ),
                              Mesh     = Angle( PhiMesh ),
                              Delta    = Angle( PhiDelta )
                              )

        if PhiOffset != 0: self.MakePhiOffset(PhiOffset)


    def MakePhiOffset(self, rotation):

        self.PhiOffset = Angle(rotation)

        r = np.ones(self.Phi.Mesh.Radian.flatten().shape)


        x, y, z = cs.sp2cart(r,
                              self.Phi.Mesh.Radian.flatten()-np.pi/2/1.01, ###########################___MAY BE A PROBLEM HERE!
                              self.Theta.Mesh.Radian.flatten()/1.01,
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
                        antialiased = False,               )

        ax.scatter(xp,yp,zp,color="k",s=2)
        ax.set_xlabel('X-Direction')
        ax.set_ylabel('Y-Direction')
        ax.set_zlabel('Z-Direction')
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])
        ax.set_zlim([-1,1])
        ax.set_aspect("auto")
        plt.tight_layout()
        plt.show()


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)




class DirectMeshes(object):

    def __init__(self,
                 Npts:     int  = 101,
                 XBound:   list = [-1,1],
                 YBound:   list = [-1,1],
                 XNpts:    int  = None,
                 YNpts:    int  = None,
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
