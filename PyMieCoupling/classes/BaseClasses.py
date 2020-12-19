import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import numpy as np

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.classes.Representations import S1S2, SPF, Stokes, Field
from PyMieCoupling.functions.Couplings import Coupling, GetFootprint
from PyMieCoupling.cpp.S1S2 import GetFieldsFromMesh
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Angle, Polarization

class BaseDetector(object):
    """Short summary.

    Parameters
    ----------
    size : float
        Size of the detector, [diameter for circle shaped/ side for square].
    shape : str
        Shape of the detector.
    wavelength : float
        Wavelength of the incoming source.
    npts : int
        Number of points defining the rastered meshes.
    ThetaOffset : float
        Offset of theta angle between the detector and source.
    PhiOffset : float
        Offset of phi angle between the detector and source.
    Magnification : float
        Magnification induced by the lense.
    Name : str
        Name of detector [optional for plots].

    Attributes
    ----------
    _name : type
        Description of attribute `_name`.
    _coupling : type
        Description of attribute `_coupling`.
    k : type
        Description of attribute `k`.
    DirectVec : type
        Description of attribute `DirectVec`.
    GenShift : type
        Description of attribute `GenShift`.
    GenMeshes : type
        Description of attribute `GenMeshes`.
    Field : type
        Description of attribute `Field`.
    Fourier : type
        Description of attribute `Fourier`.
    GenField : type
        Description of attribute `GenField`.
    magnificate : type
        Description of attribute `magnificate`.
    size
    wavelength
    ThetaOffset
    PhiOffset
    npts

    """

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
        self._Filter = Polarization(val)

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
        if val >= 1: val = 1
        if val <= 0: val = 0
        self.FarField.NA = val



    def Coupling(self, Scatterer, Mode = 'Centered'):
        return Coupling(Scatterer = Scatterer, Detector = self, Mode = Mode)


    def Footprint(self, Scatterer):
        return GetFootprint(Scatterer = Scatterer, Detector = self)



class _BaseFarField(np.ndarray):

    def __new__(cls,
                Scalar        = None,
                Meshes        = None):

        cls.Meshes = Meshes

        this = np.array(Scalar, copy=False)

        this = np.asarray(this).view(cls)

        return this


    def __array_finalize__(self, obj):
        pass


    def __init__(self, Scalar, Meshes):
        pass


    def Plot(self):

        fig, axes = plt.subplots(nrows = 1,
                                 ncols = 2,
                                 figsize    = (8,3),
                                 subplot_kw = {'projection':'mollweide'})

        axes[0].pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     self.Scalar.real,
                     shading='auto')

        axes[0].set_title('Real Part\n Far-Field Spherical Coordinates')
        axes[0].set_ylabel(r'Angle $\phi$ [Degree]')
        axes[0].set_xlabel(r'Angle $\theta$ [Degree]')

        axes[1].pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     self.Scalar.imag,
                     shading='auto')

        axes[1].set_title('Imaginary Part\n Far-Field Spherical Coordinates')
        axes[1].set_ylabel(r'Angle $\phi$ [Degree]')
        axes[1].set_xlabel(r'Angle $\theta$ [Degree]')
        axes[1].grid()

        fig.tight_layout()
        plt.show()





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
                     np.real(self.Scalar),
                     shading='auto')

        ax1.pcolormesh(
                     self.Meshes.Theta.Mesh.Radian,
                     self.Meshes.Phi.Mesh.Radian+np.pi/2,
                     np.imag(self.Scalar),
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
                                  ThetaNpts          = self.Scalar.shape[0],
                                  PhiNpts            = self.Scalar.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = 0)
    @PhiBound.setter
    def PhiBound(self, val: list):
        self._PhiBound = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Scalar.shape[0],
                                  PhiNpts            = self.Scalar.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = 0)

    @PhiOffset.setter
    def PhiOffset(self, val):
        self._PhiOffset = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Scalar.shape[0],
                                  PhiNpts            = self.Scalar.shape[1],
                                  ThetaOffset        = 0,
                                  PhiOffset          = val)

    @ThetaOffset.setter
    def ThetaOffset(self, val):
        self._ThetaOffset = val
        self.Meshes = AngleMeshes(ThetaBound         = self._ThetaBound,
                                  PhiBound           = self._PhiBound,
                                  ThetaNpts          = self.Scalar.shape[0],
                                  PhiNpts            = self.Scalar.shape[1],
                                  ThetaOffset        = val,
                                  PhiOffset          = 0)






class BaseScatterer(object):
    """Object containing all scatterer-related attributes.

    Parameters
    ----------
    diameter : float
        Diameter of the scatterer.
    wavelength : float
        Wavelength of the incident lightfield.
    index : float
        Refractive index of the scatterer.
    npts : int
        Number of points for the full solid angle of the far-field, later to
        be interpolated.

    Attributes
    ----------
    Full : <Fields class>
        It represents the entire Far-field representation of the scatterer.
    ComputeS1S2 : type
        Methode using package PyMieScatt to compute S1 and S2 parameter form mu value.
    diameter
    wavelength
    index
    npts

    """
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
        Para, Perp = GetFieldsFromMesh(m                    = self.Index,
                                       x                    = self.SizeParam,
                                       ThetaMesh            = self.Meshes.Theta.Mesh.Radian.flatten(),
                                       PhiMesh              = self.Meshes.Phi.Mesh.Radian.flatten(),
                                       Shape                = self.Meshes.Theta.Mesh.Radian.shape,
                                       Polarization         = self.Source.Polarization.Radian);

        self._FarField = Field(Perpendicular = Perp, Parallel = Para, Meshes = self.Meshes);


    def Coupling(self, Detector, Filter = None, Mode = 'Centered'):
        return Coupling(Scatterer    = self,
                        Detector     = Detector,
                        Filter       = Filter,
                        Mode         = Mode)

    def Footprint(self, Detector):
        return GetFootprint(Scatterer = self, Detector = Detector)




# -
