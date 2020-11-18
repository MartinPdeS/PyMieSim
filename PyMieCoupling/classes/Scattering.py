
import numpy as np
from PyMieCoupling.classes.Fields import Field
from PyMieCoupling.classes.Meshes import ScatMeshes
from PyMieCoupling.classes.Representations import S1S2
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.cpp.S1S2 import MieS1S2 #_CYTHON PACKAGE



class Scatterer(object):
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

    def __init__(self,
                 Diameter:    float,
                 Source:      Source,
                 Index:       float,
                 Npts:        int         = None,
                 Meshes:      ScatMeshes  = None,
                 ThetaBound:  list        = [-180, 180],
                 ThetaOffset: float       = 0,
                 PhiBound:    list        = [-180, 180],
                 PhiOffset:   float       = 0) -> None:


        self.Diameter, self.Source, self.Index = Diameter, Source, Index

        self.SizeParam = Source.k * ( self.Diameter * 2 )


        self.GenMesh(Meshes, ThetaBound, PhiBound, ThetaOffset, PhiOffset, Npts)

        self.__S1S2, self.__Field = None, None



    def GenMesh(self,
                Meshes:      ScatMeshes = None,
                ThetaBound:  list       = [-90, 90],
                PhiBound:    list       = [-90, 90],
                ThetaOffset: float      = 0,
                PhiOffset:   float      = 0,
                Npts:        int        = 101):

        if Meshes:
            self.Meshes = Meshes

        else:
            self.Meshes = ScatMeshes(ThetaBound = np.array(ThetaBound) + ThetaOffset,
                                     PhiBound   = np.array(PhiBound) + PhiOffset,
                                     Npts       = Npts)


    @property
    def S1S2(self) -> np.ndarray:
        if self.__S1S2 is None:
            self.__S1S2 = S1S2(SizeParam  = self.SizeParam,
                               Index      = self.Index,
                               Meshes     = self.Meshes)
            return self.__S1S2

        else:
            return self.__S1S2


    @property
    def Field(self) -> ScatMeshes:
        if self.__Field is None:
            self.GenField()
            return self.__Field
        else:
            return self.__Field



    def GenField(self):
        """The methode generate the <Fields> class from S1 and S2 value computed
        with the PyMieScatt package.

        """


        S1, S2 = MieS1S2(self.Index,
                         self.SizeParam,
                         self.Meshes.Phi.Vector.Radian.tolist(),
                         self.Meshes.Theta.Vector.Radian.tolist(),
                         )

        Parallel = np.outer(S1, np.sin(self.Meshes.Theta.Vector.Radian))

        Perpendicular = np.outer(S2, np.cos(self.Meshes.Theta.Vector.Radian))

        self.__Field = Field(Perpendicular = Perpendicular,
                             Parallel      = Parallel,
                             Meshes        = self.Meshes
                             )

    @property
    def Stokes(self) -> None:
        return self.Field.Stokes

    @property
    def Jones(self) -> None:
        return self.Field.Jones

    @property
    def SPF(self) -> None:
        return self.Field.SPF



    # -
