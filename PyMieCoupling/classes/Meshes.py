import numpy as np
from PyMieCoupling.functions.converts import rad2deg, deg2rad


class AngleMeshes(object):

    def __init__(self,
                 ThetaBound: list = [-180,180],
                 PhiBound:   list = [-180,180],
                 ThetaNpts:  int  = None,
                 PhiNpts:    int  = None,
                 ):

        self.ThetaNpts = ThetaNpts
        self.PhiNpts   = PhiNpts
        self.MakeProperties(ThetaBound= ThetaBound, PhiBound= PhiBound)


    def MakeProperties(self, ThetaBound, PhiBound):

        ThetaRange              = np.abs(ThetaBound[0] - ThetaBound[1])
        PhiRange                = np.abs(PhiBound[0] - PhiBound[1])

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

        self.Theta = Namespace(Boundary  = Angle( ThetaBound ),
                                Range    = Angle( ThetaRange ),
                                Vector   = Angle( ThetaVector ),
                                Mesh     = Angle( ThetaMesh ),
                                Delta    = Angle( ThetaDelta ),
                                )

        self.Phi = Namespace(Boundary  = Angle( PhiBound ),
                              Range    = Angle( PhiRange ),
                              Vector   = Angle( PhiVector ),
                              Mesh     = Angle( PhiMesh ),
                              Delta    = Angle( PhiDelta )
                              )



class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

class Angle(object):

    def __init__(self, input):
        self.Degree = input
        self.Radian = deg2rad(input)





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

        self.X = Namespace(Boundary =  XBound ,
                           Range    =  XRange ,
                           Vector   =  XVector ,
                           Mesh     =  XMesh ,
                           Delta    =  XDelta ,
                            )

        self.Y = Namespace(Boundary  =  YBound ,
                           Range    =  YRange ,
                           Vector   =  YVector ,
                           Mesh     =  YMesh ,
                           Delta    =  YDelta
                              )




# -
