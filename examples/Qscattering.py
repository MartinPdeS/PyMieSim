matplotlib=True
mlab=False

def run():
    import numpy as np
    from PyMieSim.Scatterer  import Sphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    scatKwargs   = { 'Diameter' : np.linspace(100e-9, 10000e-9, 400),
                     'Index'    : np.linspace(1.5, 1.8, 3).round(1),
                     'nMedium'  : [1] }

    sourceKwargs = { 'Wavelength'   : [400e-9],
                     'Polarization' : [0]}

    scatSet    = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

    sourceSet  = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = None)

    Eff = Experiment.Efficiencies(['Qsca'], AsType='pymiesim')

    Eff.Plot(x='Diameter')


if __name__ == '__main__':
    run()
