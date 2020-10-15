from numpy import linspace, meshgrid, pi, array, mod

from src.functions.converts import rad2deg, deg2rad


class Meshes(object):

    def __init__(self, npts, ThetaBound=[0,0], PhiBound=[0,360]):

        self.npts = npts

        ThetaBound, PhiBound = array(ThetaBound), array(PhiBound)

        self.ThetaBound = Angle(ThetaBound)

        self.PhiBound = Angle(PhiBound)

        self.MakeVec()

        self.MakeMeshes()


    def MakeVec(self):

        self.ThetaVec = Angle( linspace(*self.ThetaBound.Degree, self.npts) )

        self.PhiVec = Angle( linspace(*self.PhiBound.Degree, self.npts) )

        #self.ThetaRadVec, self.PhiRadVec = deg2rad(self.ThetaAngleVec), deg2rad(self.PhiAngleVec)


    def MakeMeshes(self):

        self.ThetaMesh, self.PhiMesh = meshgrid(self.ThetaVec.Degree, self.PhiVec.Degree)

        self.ThetaMesh = Angle(self.ThetaMesh)

        self.PhiMesh = Angle(self.PhiMesh)

        #self.ThetaAngleMesh, self.PhiAngleMesh = rad2deg(self.ThetaRadMesh), rad2deg(self.PhiRadMesh)





class Angle(object):

    def __init__(self, input):

        input = array(input)

        self.Degree = input

        self.Radian = deg2rad(input)

        self.DegreeMod = mod(self.Degree,360)

        self.RadianMod = mod(self.Radian, 2*pi)

    @property
    def Degree(self):
        return self._Degree

    @Degree.setter
    def Degree(self, Degree):
        self._Degree = Degree








# -
