from numpy import linspace, meshgrid, pi, array, mod
import numpy as np
from PyMieCoupling.functions.converts import rad2deg, deg2rad


class Meshes(object):

    def __init__(self,
                 Npts: int = 101,
                 ThetaBound: list = [-180,180],
                 PhiBound: list = [-180,180]):

        self.Npts = Npts

        ThetaBound, PhiBound = array(ThetaBound), array(PhiBound)

        self.Theta, self.Phi = _Angle(), _Angle()

        self.Theta.Boundary, self.Phi.Boundary = Angle(ThetaBound), Angle(PhiBound)

        self.MakeVec()

        self.MakeMeshes()

        self.MakeDelta()


    def MakeVec(self):
        self.Theta.Vector = Angle( linspace(*self.Theta.Boundary.Degree, self.Npts) )

        self.Phi.Vector = Angle( linspace(*self.Phi.Boundary.Degree, self.Npts) )


    def MakeMeshes(self):
        ThetaMesh, PhiMesh = meshgrid(self.Theta.Vector.Degree, self.Phi.Vector.Degree)

        self.Theta.Mesh, self.Phi.Mesh = Angle(ThetaMesh), Angle(PhiMesh)


    def MakeDelta(self):
        ThetaDelta = np.abs(self.Theta.Boundary.Degree[0] - self.Theta.Boundary.Degree[1]) / self.Npts

        PhiDelta = np.abs(self.Phi.Boundary.Degree[0] - self.Phi.Boundary.Degree[1]) / self.Npts

        self.Theta.Delta, self.Phi.Delta = Angle(ThetaDelta), Angle(PhiDelta)



class _Angle(object):

    def __init__(self):
        self.Boundary = None
        self.Vector = None
        self.Mesh = None
        self.Delta = None


class Angle(object):

    def __init__(self, input):
        self.Degree = np.array(input)

        self.Radian = deg2rad(input)

        self.DegreeMod = np.mod(self.Degree,180)

        self.RadianMod = np.mod(self.Radian, pi)


    @property
    def Degree(self):
        return self._Degree

    @Degree.setter
    def Degree(self, Degree):
        Degree = np.array(Degree)
        if (Degree < -180.1).any():
            raise Exception( 'Angle must defined between -180 to 180' )
        self._Degree = Degree








# -
