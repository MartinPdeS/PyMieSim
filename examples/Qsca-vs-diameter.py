matplotlib=True
mlab=False

def run():
    import numpy as np
    from PyMieSim            import Material
    from PyMieSim.Scatterer  import Sphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    Detector0 = Photodiode(NA                = 0.1,
                           Sampling          = 300,
                           GammaOffset       = 20,
                           PhiOffset         = 30,
                           CouplingMode      = 'Centered')

    scatKwargs   = { 'Diameter' : np.geomspace(6.36e-09, 100000e-9, 1500),
                     'Index'    : [1.4,1.8],
                     'nMedium'  : [1] }

    sourceKwargs = { 'Wavelength'   : [400e-9],
                     'Polarization' : [0]}

    scatSet   = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

    sourceSet = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    DetecSet  = DetectorSet([Detector0])

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = DetecSet)

    Qsca = Experiment.Efficiencies(Eff='Qsca', AsType='pymiesim')

    Qsca.Plot(x='Diameter', Scale='log')



if __name__ == '__main__':
    run()
