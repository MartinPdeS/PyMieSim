import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
import numpy as np
import matplotlib
import cartopy.crs as ccrs

from PyMieCoupling.classes.Meshes import AngleMeshes
from PyMieCoupling.classes.Representations import S1S2, SPF, Stokes, Field, ScalarFarField
from PyMieCoupling.functions.Couplings import Coupling, GetFootprint
from PyMieCoupling.cpp.S1S2 import GetFieldsFromMesh
from PyMieCoupling.functions.converts import NA2Angle
from PyMieCoupling.utils import Angle, _Polarization, PlotFarField, InterpFull, PlotUnstructuredSphere, PlotStructuredSphere

try:
    from PyMieCoupling.cpp.S1S2 import GetS1S2
except:
    try:
        from PyMieCoupling.cython.S1S2 import GetS1S2
    except:
        try:
            from PyMieCoupling.cython.S1S2 import GetS1S2
        except: ImportError






class MeshProperty(object):
    """Short summary.

    """

    def __init__(self):
        pass

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
        return self.Meshes.GammaOffset

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
        self.MaxAngle = NA2Angle(val).Radian
        self.Meshes.UpdateSphere(MaxAngle = self.MaxAngle)
        self.GetSpherical()




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

    """


    def Coupling(self, Scatterer):
        return Coupling(Scatterer = Scatterer, Detector = self)


    def Footprint(self, Scatterer):
        return GetFootprint(Scatterer = Scatterer, Detector = self)


    def Plot(self, num=400, scatter=True):

        return PlotUnstructuredSphere(Phi     = self.Meshes.Phi.Degree,
                                     Theta   = self.Meshes.Theta.Degree,
                                     Scalar  = self.Scalar)






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


    def S1S2(self, Num=200):
        if self._S1S2 is None:
            self._S1S2 = S1S2(Parent=self, Num=Num)
            return self._S1S2

        else:
            return self._S1S2



    def Field(self, Num=200):

        self._Field = ScalarFarField(Num = Num, Parent = self)

        return self._Field


    def Parallel(self, Phi, Theta):
        if not isinstance(self._Parallel, np.ndarray):
            self._Parallel, self._Perpendicular = self.GenField(Phi, Theta)
            return self._Parallel
        else:
            return self._Parallel


    def Perpendicular(self, Phi, Theta):
        if not isinstance(self._Perpendicular, np.ndarray):
            self._Parallel, self._Perpendicular = self.GenField(Phi, Theta)
            return self._Perpendicular
        else:
            return self._Perpendicular


    def SPF(self, Num=100):
        if not self._SPF:
            self._SPF = SPF(Parent=self, Num=Num)
            return self._SPF
        else:
            return self._SPF


    def GenField(self, Phi, Theta):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.
        """

        return GetFieldsFromMesh(m            = self.Index,
                                 x            = self.SizeParam,
                                 ThetaMesh    = Theta,
                                 PhiMesh      = Phi - np.pi/2,
                                 Polarization = self.Source.Polarization.Radian)


    def Plot(self, Num=200, scatter=False):

        Theta, Phi = np.mgrid[0:2*np.pi:complex(Num), -np.pi/2:np.pi/2:complex(Num)]

        Para, Perp = self.GenField(Phi.flatten(), Theta.flatten())

        fig0 = PlotStructuredSphere(Phi     = np.rad2deg(Phi),
                                    Theta   = np.rad2deg(Theta),
                                    Scalar  = Para.reshape(Theta.shape))

        fig1 = PlotStructuredSphere(Phi     = np.rad2deg(Phi),
                                    Theta   = np.rad2deg(Theta),
                                    Scalar  = Perp.reshape(Theta.shape))

        return fig0, fig1


    def Footprint(self, Detector):
        return GetFootprint(Scatterer = self, Detector = Detector)






# -
