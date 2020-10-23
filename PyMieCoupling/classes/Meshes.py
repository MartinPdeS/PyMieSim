from numpy import linspace, meshgrid, pi, array, mod
import numpy as np
from PyMieCoupling.functions.converts import rad2deg, deg2rad


class Meshes(object):

    def __init__(self,
                 npts: int = 101,
                 ThetaBound: list = [-180,180],
                 PhiBound: list = [-180,180]):

        self.npts = npts

        ThetaBound, PhiBound = array(ThetaBound), array(PhiBound)

        self.Theta, self.Phi = _Angle(), _Angle()

        self.Theta.Boundary, self.Phi.Boundary = Angle(ThetaBound), Angle(PhiBound)

        self.MakeVec()

        self.MakeMeshes()


    def MakeVec(self):
        self.Theta.Vector = Angle( linspace(*self.Theta.Boundary.Degree, self.npts) )

        self.Phi.Vector = Angle( linspace(*self.Phi.Boundary.Degree, self.npts) )


    def MakeMeshes(self):
        ThetaMesh, PhiMesh = meshgrid(self.Theta.Vector.Degree, self.Phi.Vector.Degree)

        self.Theta.Mesh, self.Phi.Mesh = Angle(ThetaMesh), Angle(PhiMesh)


class _Angle(object):

    def __init__(self):
        self.Boundary = None
        self.Vector = None
        self.Mesh = None


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
