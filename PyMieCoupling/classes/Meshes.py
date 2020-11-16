import numpy as np
from PyMieCoupling.functions.converts import rad2deg, deg2rad
from PyMieCoupling.classes.Misc import Operation as Op

class Meshes(object):

    def __init__(self,
                 Npts:       int  = 101,
                 ThetaBound: list = [-180,180],
                 PhiBound:   list = [-180,180]):

        self.Npts = Npts

        self.MakeBound(ThetaBound= ThetaBound, PhiBound= PhiBound)

        self.MakeVec()

        self.MakeMeshes()

        self.MakeDeltas()


    def MakeBound(self, ThetaBound, PhiBound):
        self.Theta, self.Phi = _Angle(), _Angle()

        self.Theta.Boundary = Angle( ThetaBound )

        self.Phi.Boundary = Angle( PhiBound )

        self.Theta.Range = Angle( np.abs(ThetaBound[0] - ThetaBound[1]) )

        self.Phi.Range = Angle( np.abs(PhiBound[0] - PhiBound[1]) )


    def MakeVec(self):
        self.Theta.Vector = Angle( np.linspace(*self.Theta.Boundary.Degree, self.Npts) )

        self.Phi.Vector = Angle( np.linspace(*self.Phi.Boundary.Degree, self.Npts) )


    def MakeMeshes(self):
        ThetaMesh, PhiMesh = np.meshgrid(self.Theta.Vector.Degree, self.Phi.Vector.Degree)

        self.Theta.Mesh, self.Phi.Mesh = Angle(ThetaMesh), Angle(PhiMesh)


    def MakeDeltas(self):
        ThetaDelta = (self.Theta.Boundary.Degree[0] - self.Theta.Boundary.Degree[1]).__abs__() / self.Npts

        PhiDelta = (self.Phi.Boundary.Degree[0] - self.Phi.Boundary.Degree[1]).__abs__() / self.Npts

        self.Theta.Delta, self.Phi.Delta = Angle(ThetaDelta), Angle(PhiDelta)




class _Angle(object):

    def __init__(self):
        self.Boundary = None
        self.Vector = None
        self.Mesh = None
        self.Delta = None
        self.Range = None


class Angle(object):

    def __init__(self, input):
        self.Degree = input
        self.Radian = deg2rad(input)









# -
