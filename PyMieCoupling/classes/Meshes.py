import numpy as np
import cupy as cp
from PyMieCoupling.functions.converts import rad2deg, deg2rad


class Meshes(object):

    def __init__(self,
                 Npts:       int  = 101,
                 ThetaBound: list = [-180,180],
                 PhiBound:   list = [-180,180],
                 GPU:        bool = False):

        self.Npts = Npts

        self.GPU = GPU

        self.MakeBound(ThetaBound= ThetaBound, PhiBound= PhiBound)

        self.MakeVec()

        self.MakeMeshes()

        self.MakeDeltas()


    def MakeBound(self, **kwargs):
        self.Theta, self.Phi = _Angle(), _Angle()

        if self.GPU: return self.MakeBoundGPU(**kwargs)

        else: return self.MakeBoundCPU(**kwargs)


    def MakeVec(self):
        if self.GPU: return self.MakeVecGPU()

        else: return self.MakeVecCPU()


    def MakeMeshes(self):
        if self.GPU: return self.MakeMeshesGPU()

        else: return self.MakeMeshesCPU()


    def MakeDeltas(self):
        if self.GPU: self.MakeDeltasGPU()

        else: return self.MakeDeltasCPU()


    def MakeBoundCPU(self, ThetaBound, PhiBound):
        self.Theta.Boundary = Angle( np.array(ThetaBound) )

        self.Phi.Boundary =  Angle( np.array(PhiBound) )


    def MakeBoundGPU(self, ThetaBound, PhiBound):
        self.Theta.Boundary = Angle( cp.array(ThetaBound) )

        self.Phi.Boundary = Angle( cp.array(PhiBound) )


    def MakeVecCPU(self):
        self.Theta.Vector = Angle( np.linspace(*self.Theta.Boundary.Degree, self.Npts) )

        self.Phi.Vector = Angle( np.linspace(*self.Phi.Boundary.Degree, self.Npts) )


    def MakeVecGPU(self):
        self.Theta.Vector = Angle( cp.linspace(*self.Theta.Boundary.Degree, self.Npts) )

        self.Phi.Vector = Angle( cp.linspace(*self.Phi.Boundary.Degree, self.Npts) )


    def MakeMeshesCPU(self):
        ThetaMesh, PhiMesh = np.meshgrid(self.Theta.Vector.Degree, self.Phi.Vector.Degree)

        self.Theta.Mesh, self.Phi.Mesh = Angle(ThetaMesh), Angle(PhiMesh)


    def MakeMeshesGPU(self):
        ThetaMesh, PhiMesh = cp.meshgrid(self.Theta.Vector.Degree, self.Phi.Vector.Degree)

        self.Theta.Mesh, self.Phi.Mesh = Angle(ThetaMesh), Angle(PhiMesh)


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
