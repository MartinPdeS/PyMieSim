import numpy as np
import os.path
from typing   import Union
from beartype import beartype


import PyMieSim
from PyMieSim.BaseClasses import BaseDetector, MeshProperty
from PyMieSim.Physics     import _Polarization, Angle
from PyMieSim.utils       import ( interp_at,
                                   NA2Angle,
                                   Normalize,
                                   RescaleComplex,
                                   RotateComplex )


class Photodiode(BaseDetector, MeshProperty):
    """Detector type class representing a photodiode, light coupling is
    thus independant of the phase of the latter.

    Parameters
    ----------
    NA : :class:`float`
        Numerical aperture of imaging system.
    Sampling : :class:`int`
        Number of sampling points for the mode (inside NA).
    GammaOffset : :class:`float`
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : :class:`float`
        Angle offset of detector in the direction parallel to polarization.
    Filter : :class:`float`
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : :class:`str`
        Methode for computing mode coupling. Either Centered or Mean.

    """
    @beartype
    def __init__(self,
                 NA           : Union[int, float],
                 Sampling     : int                      = 400,
                 CouplingMode : str                      = 'Centered',
                 GammaOffset  : Union[int, float]        = 0,
                 PhiOffset    : Union[int, float]        = 0,
                 Filter       : Union[int, float, bool]  = None,
                 Testing      : bool                     = False,):

        self.Testing = Testing

        self.CouplingMode = ('Intensity', CouplingMode)

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self._NA = NA2Angle(NA).Radian

        self.Mesh = self.SphericalMesh(Sampling    = Sampling,
                                       MaxAngle    = self._NA,
                                       PhiOffset   = PhiOffset,
                                       GammaOffset = GammaOffset,
                                       Structured  = False)

        self.Scalar = self.FarField(Num=self.Mesh.Sampling, Structured=False)


    def FarField(self, Num=251, Structured=False):
        """
        Compute the FarField in a structured or unstructured Mesh.

        Parameters
        ----------
        Num : :class:`int`
            Dimension of the structured mesh [Num, Num].

        Structured : :class:`bool`
            Indicate how to compute the far field.

        Returns
        -------
        :class`numpy.ndarray`:
            FarField array.

        """
        if Structured: return np.ones([Num, Num])
        else:          return np.ones(Num)


    def __repr__(self):

        return f"""
        Photodiode detector
        Coupling Mode: Intensity
        Sampling:      {self.Mesh.Sampling}
        Gamma offset:  {self.Mesh.GammaOffset}
        Phi offset:    {self.Mesh.PhiOffset}
        """




class IntegratingSphere(Photodiode):
    """Detector type class representing a photodiode, light coupling is
    thus independant of the phase of the latter.

    Parameters
    ----------
    NA : :class:`float`
        Numerical aperture of imaging system.
    Sampling : :class:`int`
        Number of sampling points for the mode (inside NA).
    GammaOffset : :class:`float`
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : :class:`float`
        Angle offset of detector in the direction parallel to polarization.
    Filter : :class:`float`
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : :class:`str`
        Methode for computing mode coupling. Either Centered or Mean.

    """
    @beartype
    def __init__(self,
                 Sampling     : int                      = 400,
                 CouplingMode : str                      = 'Centered',
                 Filter       : Union[int, float, bool]  = None,
                 Testing      : bool                     = False,):

        self.Testing = Testing

        self.CouplingMode = ('Intensity', CouplingMode)

        self._GammaOffset, self._PhiOffset = 0, 0

        self._Filter = _Polarization(Filter)

        self._NA = 2.0

        self.Mesh = self.SphericalMesh(Sampling    = Sampling,
                                       MaxAngle    = self._NA,
                                       PhiOffset   = 0,
                                       GammaOffset = 0,
                                       Structured  = False)

        self.Scalar = self.FarField(Num=self.Mesh.Sampling, Structured=False)


    def __repr__(self):

        return f"""
        Integrating sphere
        Coupling Mode: Intensity
        Sampling:      {self.Mesh.Sampling}
        """




class LPmode(BaseDetector, MeshProperty):
    """Detector type class representing a fiber LP mode, light coupling is
    thus dependant of the phase of the latter.

    Parameters
    ----------
    Mode : :class:`tuple`
        LP mode index l, m.
    NA : :class:`float`
        Numerical aperture of imaging system.
    Sampling : :class:`int`
        Number of sampling points for the mode (inside NA).
    InterpSampling : :class:`int`
        Number of sampling point for interpolation of FarField mode.
    GammaOffset : :class:`float`
        Angle offset of detector in the direction perpendicular to polarization.
    PhiOffset : :class:`float`
        Angle offset of detector in the direction parallel to polarization.
    Filter : :class:`float`
        Angle of polarization filter in front of detector. Default is "None"
    CouplingMode : :class:`str`
        Methode for computing mode coupling. Either [Centered or Mean].

    """

    @beartype
    def __init__(self,
                 Mode         : Union[tuple, list],
                 NA           : float,
                 Rotation     : Union[int, float]        = 0,
                 Sampling     : int                      = 401,
                 GammaOffset  : Union[int, float]        = 0,
                 PhiOffset    : Union[int, float]        = 0,
                 Filter       : Union[int, float, bool]  =  None,
                 CouplingMode : str                      = 'Centered',
                 Testing      : bool                     = False):

        self.Testing = Testing

        assert CouplingMode in ['Centered','Mean'], "Coupling mode can either be Centered or Mean"

        if NA > 1 or NA < 0: print("WARNING: High values of NA do not comply \
                                   with paraxial approximation. Value \
                                   under 0.4 are prefered")

        self.CouplingMode = ('Amplitude', CouplingMode)

        self._Filter      = _Polarization(Filter)

        self.ModeNumber   = Mode

        self.Mesh = self.SphericalMesh(Sampling    = Sampling,
                                       MaxAngle    = NA2Angle(NA).Radian,
                                       PhiOffset   = PhiOffset,
                                       GammaOffset = GammaOffset,
                                       Structured  = False)

        self.Scalar = self.FarField(Num        = self.Mesh.Sampling,
                                    Structured = False,
                                    Rotation   = Rotation)


    def FarField(self, Num=251, Structured=False, Rotation=0):
        """
        Compute the FarField in a structured or unstructured Mesh.
        The unstructured far field is computed using linear interpolation
        on a structured mesh.

        Parameters
        ----------
        Num : :class:`int`
            Dimension of the structured mesh [Num, Num].

        Structured : :class:`bool`
            Indicate how to compute the far field.

        Returns
        -------
        :class`numpy.ndarray`:
            FarField array.

        """

        filename = PyMieSim.__path__[0] + f'/LPmodes/FLP{self.ModeNumber[0]}{self.ModeNumber[1]}.npy'

        if not os.path.exists(filename):
            raise ValueError("The LP mode has not been previously compilated. "
                              "Please consult the documentation to do so. "
                              "Doc available at: "
                              "https://pymiesim.readthedocs.io/en/latest/Intro.html")

        mode = np.load(filename)

        if Rotation !=0: mode = RotateComplex(mode, Rotation)

        if Num != mode.shape[0]: mode = RescaleComplex(Input=mode, Scale=Num/mode.shape[0])

        if Structured: return mode

        self._phi, self._theta = self.SphericalMesh(Sampling    = mode.shape[0],
                                                    MaxAngle    = self.Mesh.MaxAngle,
                                                    PhiOffset   = 0,
                                                    GammaOffset = 0,
                                                    Structured  = True)

        Interp = interp_at(x           = self._phi.flatten(),
                           y           = self._theta.flatten(),
                           v           = mode.astype(np.complex).flatten(),
                           xp          = self.Mesh.base[0],
                           yp          = self.Mesh.base[1],
                           algorithm   = 'linear',
                           extrapolate = True)

        return Normalize(Interp)





    def __repr__(self):
        return f"""
        LP mode detector
        Coupling Mode: Amplitude
        LP Mode:       {(self.ModeNumber[0]-1, self.ModeNumber[1], self.ModeNumber[2])}
        Sampling:      {self.Mesh.Sampling}
        Gamma offset:  {self.Mesh.GammaOffset}
        Phi offset:    {self.Mesh.PhiOffset}
        """






class _Photodiode(BaseDetector, MeshProperty):
    """Detector class for develop use only. Do not use!

    """

    @beartype
    def __init__(self,
                 NA           : Union[int, float],
                 Sampling     : int                      = 400,
                 CouplingMode : str                      = 'Centered',
                 GammaOffset  : Union[int, float]        = 0,
                 PhiOffset    : Union[int, float]        = 0,
                 Filter       : Union[int, float, bool]  = None):

        self.CouplingMode = ('Intensity', CouplingMode)

        self._GammaOffset, self._PhiOffset = GammaOffset, PhiOffset

        self._Filter = _Polarization(Filter)

        self.Mesh = self.SphericalMesh(Sampling    = Sampling,
                                       MaxAngle    = NA2Angle(NA).Radian,
                                       PhiOffset   = PhiOffset,
                                       GammaOffset = GammaOffset,
                                       Structured  = False)

        self.Scalar = self.FarField(Num=self.Mesh.Sampling, Structured=False)


    def FarField(self, Num=251, Structured=False):
        """
        Compute the FarField in a structured or unstructured Mesh.

        Parameters
        ----------
        Num : :class:`int`
            Dimension of the structured mesh [Num, Num].

        Structured : :class:`bool`
            Indicate how to compute the far field.

        Returns
        -------
        :class`numpy.ndarray`:
            FarField array.

        """
        if Structured: return np.ones([Num, Num])
        else:          return np.ones(Num)


    def __repr__(self):

        return f"""
        Photodiode detector
        Coupling Mode: Intensity
        Sampling:      {self.Mesh.Sampling}
        Gamma offset:  {self.Mesh.GammaOffset}
        Phi offset:    {self.Mesh.PhiOffset}
        """



# -
