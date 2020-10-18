from numpy import linspace, meshgrid, pi, array, mod
import numpy as np
from miecoupling.src.functions.converts import rad2deg, deg2rad


class Meshes(object):

    def __init__(self,
                 npts: int = 101,
                 ThetaBound: list = [0,0],
                 PhiBound: list = [0,360]):

        self.npts = npts

        ThetaBound, PhiBound = array(ThetaBound), array(PhiBound)

        self.ThetaBound = Angle(ThetaBound)

        self.PhiBound = Angle(PhiBound)

        self.MakeVec()

        self.MakeMeshes()


    def MakeVec(self):

        self.ThetaVec = Angle( linspace(*self.ThetaBound.Degree, self.npts) )

        self.PhiVec = Angle( linspace(*self.PhiBound.Degree, self.npts) )


    def MakeMeshes(self):

        self.ThetaMesh, self.PhiMesh = meshgrid(self.ThetaVec.Degree, self.PhiVec.Degree)

        self.ThetaMesh = Angle(self.ThetaMesh)

        self.PhiMesh = Angle(self.PhiMesh)


class Angle(object):

    def __init__(self, input):

        input = array(input)

        self.Degree = input

        self.Radian = deg2rad(input)

        self.DegreeMod = mod(self.Degree,180)

        self.RadianMod = mod(self.Radian, pi)


    @property
    def Degree(self):
        return self._Degree

    @Degree.setter
    def Degree(self, Degree):
        Degree = np.array(Degree)
        if (Degree < -180).any():
            raise Exception( 'Angle must defined between -180 to 180' )
        self._Degree = Degree








# -
