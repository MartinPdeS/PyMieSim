matplotlib=True
mlab=False

def run():
    import numpy as np
    from PyMieSim.Scatterer  import Sphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    Detector0 = Photodiode(NA                = 0.1,
                           Sampling          = 300,
                           GammaOffset       = 20,
                           PhiOffset         = 30,
                           CouplingMode      = 'Centered')

    scatKwargs   = { 'Diameter' : np.linspace(400e-9, 1000e-9, 100),
                     'Index'    : [1.5],
                     'nMedium'  : [1,1.3] }

    sourceKwargs = { 'Wavelength'   : [400e-9],
                     'Polarization' : [0]}

    scatSet   = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

    sourceSet = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    DetecSet  = DetectorSet([Detector0])

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = DetecSet)

    Qsca = Experiment.Coupling(AsType='pymiesim')

    Qsca.Plot(x='Diameter')


if __name__ == '__main__':
    run()
