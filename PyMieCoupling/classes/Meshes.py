import numpy as np
from PyMieCoupling.functions.converts import rad2deg, deg2rad

class ScatMeshes(object):

    def __init__(self,
                 Npts:       int  = 101,
                 ThetaBound: list = [-180,180],
                 PhiBound:   list = [-180,180]):

        self.Npts = Npts

        self.MakeProperties(ThetaBound= ThetaBound, PhiBound= PhiBound)


    def MakeProperties(self, ThetaBound, PhiBound):

        ThetaRange              = np.abs(ThetaBound[0] - ThetaBound[1])
        PhiRange                = np.abs(PhiBound[0] - PhiBound[1])

        ThetaVector             = np.linspace(ThetaBound[0], ThetaBound[1], self.Npts)
        PhiVector               = np.linspace(*PhiBound, self.Npts)

        ThetaMesh, PhiMesh      = np.meshgrid(ThetaVector, PhiVector)

        ThetaDelta = np.abs(ThetaBound[0] - ThetaBound[1]) / self.Npts

        PhiDelta = np.abs(PhiBound[0] - PhiBound[1]) / self.Npts

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









# -
