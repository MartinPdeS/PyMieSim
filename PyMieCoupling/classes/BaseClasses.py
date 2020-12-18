
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



# -
