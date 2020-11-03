import numpy as np
import cupy as cp
from PyMieCoupling.functions.converts import rad2deg, deg2rad
from PyMieCoupling.classes.Misc import Operation as Op

class Meshes(object):

    def __init__(self,
                 Npts:       int  = 101,
                 ThetaBound: list = [-180,180],
                 PhiBound:   list = [-180,180],
                 cuda:       bool = False):

        self.Npts = Npts

        self.cuda = cuda

        self.MakeBound(ThetaBound= ThetaBound, PhiBound= PhiBound)

        self.MakeVec()

        self.MakeMeshes()

        self.MakeDeltas()


    def MakeBound(self, ThetaBound, PhiBound):
        self.Theta, self.Phi = _Angle(), _Angle()

        self.Theta.Boundary = Angle( ThetaBound )

        self.Phi.Boundary = Angle( PhiBound )


    def MakeVec(self):
        self.Theta.Vector = Angle( Op.linspace(self.cuda)(*self.Theta.Boundary.Degree, self.Npts) )

        self.Phi.Vector = Angle( Op.linspace(self.cuda)(*self.Phi.Boundary.Degree, self.Npts) )


    def MakeMeshes(self):
        ThetaMesh, PhiMesh = Op.meshgrid(self.cuda)(self.Theta.Vector.Degree, self.Phi.Vector.Degree)

        self.Theta.Mesh, self.Phi.Mesh = Angle(ThetaMesh), Angle(PhiMesh)


    def MakeDeltas(self):
        if self.cuda: self.MakeDeltasGPU()

        else: return self.MakeDeltasCPU()


    def MakeDeltasCPU(self):
        ThetaDelta = np.abs(self.Theta.Boundary.Degree[0] - self.Theta.Boundary.Degree[1]) / self.Npts

        PhiDelta = np.abs(self.Phi.Boundary.Degree[0] - self.Phi.Boundary.Degree[1]) / self.Npts

        self.Theta.Delta, self.Phi.Delta = Angle(ThetaDelta), Angle(PhiDelta)


    def MakeDeltasGPU(self):
        ThetaDelta = cp.abs(self.Theta.Boundary.Degree[0] - self.Theta.Boundary.Degree[1]) / self.Npts

        PhiDelta = cp.abs(self.Phi.Boundary.Degree[0] - self.Phi.Boundary.Degree[1]) / self.Npts

        self.Theta.Delta, self.Phi.Delta = Angle(ThetaDelta), Angle(PhiDelta)



class _Angle(object):

    def __init__(self):
        self.Boundary = None
        self.Vector = None
        self.Mesh = None
        self.Delta = None


class Angle(object):

    def __init__(self, input):
        self.Degree = input
        self.Radian = deg2rad(input)









# -
