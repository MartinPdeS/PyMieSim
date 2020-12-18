import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import numpy as np

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.classes.Representations import S1S2, SPF, Stokes, Field
from PyMieCoupling.cpp.S1S2 import GetFields as Fields_CPP

class BaseDetector(object):
    def __init__(self):
        pass


    @property
    def ThetaBound(self):
        return self.__ThetaBound

    @ThetaBound.setter
    def ThetaBound(self, val: list):
        self.FarField.ThetaBound = val

    @property
    def Filter(self):
        return self._Filter

    @Filter.setter
    def Filter(self, val):
        self._Filter = Angle(val)

    @property
    def PhiBound(self):
        return self.FarField.__PhiBound

    @PhiBound.setter
    def PhiBound(self, val: list):
        self.FarField.PhiBound = val

    @property
    def PhiOffset(self):
        return self.FarField.__PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.FarField.PhiOffset = val

    @property
    def ThetaOffset(self):
        return self.FarField.__ThetaOffset

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self.FarField.ThetaOffset = val

    @property
    def NA(self):
        return self.FarField._NA

    @NA.setter
    def NA(self, val):
        if val >= 1:
            val = 1
        if val <= 0:
            val = 0
        self.FarField.NA = val



    def Coupling(self,
                 Scatterer,
                 Mode         = 'Centered'):

        return Coupling(Scatterer    = Scatterer,
                        Detector     = self,
                        Mode         = Mode)


    def Footprint(self, Scatterer):
        return GetFootprint(Scatterer    = Scatterer,
                            Detector     = self)






class BaseFarField(object):
    def __init__(self):
        pass


    def Plot(self):
        fig = plt.figure(figsize=(12,3))
        ax0 = fig.add_subplot(121, projection = 'mollweide')
        ax1 = fig.add_subplot(122, projection = 'mollweide')

        ax0.pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.real(self.Spherical).T,
                     shading='auto')

        ax1.pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.imag(self.Spherical).T,
                     shading='auto')

        ax0.set_title('Real Part\n Far-Field Spherical Coordinates')
        ax0.set_ylabel(r'Angle $\phi$ [Degree]')
        ax0.set_xlabel(r'Angle $\theta$ [Degree]')
        ax0.grid()

        ax1.set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        ax1.set_ylabel(r'Angle $\phi$ [Degree]')
        ax1.set_xlabel(r'Angle $\theta$ [Degree]')
        ax1.grid()

        fig.tight_layout()


    @property
    def ThetaBound(self):
        return self._ThetaBound

    @property
    def PhiBound(self):
        return self._PhiBound

    @property
    def PhiOffset(self):
        return self.__PhiOffset

    @property
    def ThetaOffset(self):
        return self.__ThetaOffset

    @property
    def NA(self):
        return self._NA

    @NA.setter
    def NA(self, val: float):
        self._NA = val
        self.PhiBound =  np.asarray( [0, NA2Angle(self._NA)] )


    @ThetaBound.setter
    def ThetaBound(self, val: list):
        self._ThetaBound = np.asarray( val )
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = 0)
    @PhiBound.setter
    def PhiBound(self, val: list):
        self._PhiBound = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = 0)

    @PhiOffset.setter
    def PhiOffset(self, val):
        self._PhiOffset = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = val)

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self._ThetaOffset = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Spherical.shape[0],
                                  PhiNpts            = self.Spherical.shape[1],
                                  ThetaOffset        = val,
                                  PhiOffset          = 0)






class BaseScatterer(object):
    def __init__(self):
        pass


    @property
    def FarField(self) -> AngleMeshes:
        if self._FarField is None:
            self.GenField()
            return self._FarField
        else:

            return self._FarField


    @property
    def S1S2(self) -> np.ndarray:
        if self._S1S2 is None:
            self._S1S2 = S1S2(SizeParam  = self.SizeParam,
                                  Index      = self.Index,
                                  Meshes     = self.Meshes)
            return self._S1S2

        else:
            return self._S1S2


    @property
    def Stokes(self) -> None:
        if not self._Stokes:
            self._Stokes = Stokes(Field = self.Field)
            return self._Stokes
        else:
            return self._Stokes


    @property
    def SPF(self) -> None:
        if not self._SPF:
            self._SPF = SPF(Index=self.Index, SizeParam=self.SizeParam)
            return self._SPF
        else:
            return self._SPF





    def GenField(self):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.
        """

        Parallel, Perpendicular = Fields_CPP(self.Index,
                                             self.SizeParam,
                                             self.Meshes.Theta.Mesh.Radian.flatten(),
                                             self.Meshes.Phi.Mesh.Radian.flatten(),
                                             self.Meshes.Phi.Vector.Radian,
                                             self.Meshes.Theta.Mesh.Radian.shape[0],
                                             self.Meshes.Theta.Mesh.Radian.shape[1],
                                             Polarization  = self.Source.Polarization.Radian);


        self._FarField = Field(Perpendicular = Perpendicular,
                                Parallel      = Parallel,
                                Meshes        = self.Meshes);





    def Coupling(self,
                 Detector,
                 Filter       = None,
                 Mode         = 'Centered'):

        return Coupling(Scatterer    = self,
                        Detector     = Detector,
                        Filter       = Filter,
                        Mode         = Mode)

    def Footprint(self, Detector):

        return GetFootprint(Scatterer    = self,
                            Detector     = Detector)




# -
