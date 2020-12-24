import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import numpy as np
import weakref

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.classes.Representations import S1S2, SPF, Stokes, Field, ScalarFarField
from PyMieCoupling.functions.Couplings import Coupling, GetFootprint
from PyMieCoupling.cpp.S1S2 import GetFieldsFromMesh
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Angle, _Polarization

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
    GammaOffset : float
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
    GammaOffset
    PhiOffset
    npts

    """

    def __init__(self):
        pass

    def Plot(self, *args, **kwargs):
        return self.Scalar.Plot(*args, **kwargs)

    @property
    def Filter(self):
        return self._Filter

    @Filter.setter
    def Filter(self, val):
        self._Filter = _Polarization(val)

    @property
    def PhiOffset(self):
        return self.Meshes.PhiOffset

    @PhiOffset.setter
    def PhiOffset(self, val):
        self.Meshes.UpdateSphere(PhiOffset = val)
        self.GetSpherical()


    @property
    def GammaOffset(self):
        return self.FarField.__GammaOffset

    @GammaOffset.setter
    def GammaOffset(self, val):
        self.Meshes.UpdateSphere(GammaOffset = val)
        self.GetSpherical()

    @property
    def NA(self):
        return self._NA

    @NA.setter
    def NA(self, val):
        if val >= 1: val = 1
        if val <= 0: val = 0
        MaxAngle = NA2Angle(val).Radian
        self.Meshes.UpdateSphere(MaxAngle = MaxAngle)
        self.GetSpherical()


    def Coupling(self, Scatterer, Mode = 'Centered'):
        return Coupling(Scatterer = Scatterer, Detector = self, Mode = Mode)


    def Footprint(self, Scatterer):
        return GetFootprint(Scatterer = Scatterer, Detector = self)







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
    def S1S2(self) -> np.ndarray:
        if self._S1S2 is None:
            self._S1S2 = S1S2(SizeParam  = self.SizeParam,
                              Index      = self.Index,
                              Meshes     = self.Meshes)
            return self._S1S2

        else:
            return self._S1S2


    @property
    def Parallel(self):
        if not isinstance(self._Parallel, np.ndarray):
            self.GenField()
            return self._Parallel
        else:
            return self._Parallel

    @property
    def Perpendicular(self):
        if not isinstance(self._Perpendicular, np.ndarray):
            self.GenField()
            return self._Perpendicular
        else:
            return self._Perpendicular


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
            self._SPF = SPF(Parent=self)
            return self._SPF
        else:
            return self._SPF



    def GenField(self):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.
        """

        self._Parallel, self._Perpendicular = GetFieldsFromMesh(m            = self.Index,
                                                                x            = self.SizeParam,
                                                                ThetaMesh    = self.Meshes.Theta.Radian,
                                                                PhiMesh      = self.Meshes.Phi.Radian,
                                                                Polarization = self.Source.Polarization.Radian);



    def Plot(self, num=200, scatter=True):
        ThetaMesh, PhiMesh = np.mgrid[-np.pi:np.pi:complex(num),
                                      -np.pi/2:np.pi/2:complex(num)]

        Para, Perp = GetFieldsFromMesh(m                    = self.Index,
                                       x                    = self.SizeParam,
                                       ThetaMesh            = ThetaMesh.flatten(),
                                       PhiMesh              = PhiMesh.flatten(),
                                       Polarization         = 0);

        fig, axes = plt.subplots(nrows    = 2,
                               ncols      = 2,
                               figsize    = (8, 6),
                               subplot_kw = {'projection':'mollweide'}
                               )
        dic = {'Parallel': Para, 'Perpendicular': Perp}
        n = 0; axes = axes.ravel()

        for key, value in dic.items():
            im = axes[n].pcolormesh(ThetaMesh,
                               PhiMesh,
                               value.real.reshape([num,num]),
                               cmap = 'RdBu',
                               shading='auto')

            plt.colorbar(mappable=im, orientation='horizontal', ax=axes[n])
            axes[n].grid()
            axes[n].set_title('Real Part [{0}] Field'.format(key))
            axes[n].set_ylabel(r'Angle $\phi$ [Degree]')
            axes[n].set_xlabel(r'Angle $\theta$ [Degree]')

            n += 1

            im = axes[n].pcolormesh(ThetaMesh,
                               PhiMesh,
                               value.imag.reshape([num,num]),
                               cmap = 'RdBu',
                               shading='auto')

            plt.colorbar(mappable=im, orientation='horizontal', ax=axes[n])
            axes[n].grid()
            axes[n].set_title('Imaginary Part [{0}] Field'.format(key))
            axes[n].set_ylabel(r'Angle $\phi$ [Degree]')
            axes[n].set_xlabel(r'Angle $\theta$ [Degree]')
            n +=1


        if scatter:
            for i in range(4):
                axes[i].scatter(self.Meshes.Theta.Radian, self.Meshes.Phi.Radian, s=0.1,c='k')


        fig.tight_layout()



    def Coupling(self, Detector, Filter = None, Mode = 'Centered'):
        return Coupling(Scatterer    = self,
                        Detector     = Detector,
                        Filter       = Filter,
                        Mode         = Mode)

    def Footprint(self, Detector):
        return GetFootprint(Scatterer = self, Detector = Detector)




# -
