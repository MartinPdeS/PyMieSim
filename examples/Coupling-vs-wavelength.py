def run():
    import numpy as np
    from PyMieSim            import Material
    from PyMieSim.Scatterer  import Sphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    scatKwargs   = { 'Diameter'    : 200e-9,
                     'Material'    : Material('BK7'),
                     'nMedium'     : [1] }

    sourceKwargs = { 'Wavelength'   : np.linspace(400e-9, 1000e-9, 500),
                     'Polarization' : [0]}

    Detector0 = Photodiode(NA                = 2.0,
                           Sampling          = 300,
                           GammaOffset       = 0,
                           PhiOffset         = 0,
                           CouplingMode      = 'Centered')

    detecSet   = DetectorSet([Detector0])

    scatSet    = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

    sourceSet  = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = detecSet)

    Coupling = Experiment.Get('Coupling')
    #print(Experiment.config)

    Coupling.Plot(x='wavelength')


if __name__ == '__main__':
    run()
