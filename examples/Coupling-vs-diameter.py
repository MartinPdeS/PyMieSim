def run():
    import numpy as np
    from PyMieSim            import Material
    from PyMieSim.Scatterer  import Sphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    scatKwargs   = { 'Diameter'    : np.linspace(400e-9, 2000e-9, 200),
                     'Material'    : Material('BK7'),
                     'nMedium'     : [1] }

    sourceKwargs = { 'Wavelength'   : 1e-6,
                     'Polarization' : [0,30,60,90]}

    Detector0 = Photodiode(NA                = 0.2,
                           Sampling          = 300,
                           GammaOffset       = 70,
                           PhiOffset         = 0,
                           CouplingMode      = 'Centered')

    detectKwargs = { 'Detector 0'   : Detector0}

    detecSet   = DetectorSet(kwargs = detectKwargs)

    scatSet    = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

    sourceSet  = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = detecSet)

    Coupling = Experiment.Get('Coupling')

    Coupling.Plot(x='Diameter')


if __name__ == '__main__':
    run()
