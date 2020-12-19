import numpy as np
import fibermodes

class Source(object):

    def __init__(self,
                 Wavelength:   float,
                 Polarization: float,
                 Power:        float = 1):

        self.Wavelength = Wavelength

        self.k = 2 * np.pi / Wavelength

        self.Power = Power

        if Polarization != None:
            self.Polarization = Angle(Polarization)
        else:
            self.Polarization = None


class Polarization(object):

    def __init__(self, input,):
        if input == 'None':
            self.Degree = 'None'
            self.Radian = 'None'
        else:
            self.Degree = input
            self.Radian = np.deg2rad(input)



class Angle(object):

    def __init__(self, input, unit='Degree'):
        if input == 'None':
            self.Degree = 'None'
            self.Radian = 'None'

        if unit == 'Degree':
            self.Degree = input
            self.Radian = np.deg2rad(input)
        if unit == 'Radian':
            self.Degree = np.rad2deg(input)
            self.Radian = input


def SMF28():
    CoreDiameter = 8.2e-6
    cladDiameter = 125e-6

    Fiber = fiber(core_radius = CoreDiameter,
                  core_index  = 1.4456,
                  clad_radius = cladDiameter,
                  clad_index  = 1.4444)

    return Fiber, CoreDiameter



class fiber(object):

    def __init__(self,
                 core_radius,
                 core_index,
                 clad_radius,
                 clad_index):

        self.MaxDirect = 2 * clad_radius

        factory = fibermodes.FiberFactory()

        factory.addLayer(name     = 'core',
                         radius   = core_radius,
                         material = 'Fixed',
                         geometry = "StepIndex",
                         index    = 1.4489)

        factory.addLayer(name     = 'cladding',
                         material = 'Fixed',
                         index    = 1)

        self.source = factory[0]
