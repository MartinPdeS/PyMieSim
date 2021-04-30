matplotlib=True
mlab=False

def run():
    import numpy as np
    from PyMieSim.Scatterer  import Sphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    scatKwargs   = { 'Diameter' : 200e-9,
                     'Index'    : [4],
                     'nMedium'  : [1] }

    sourceKwargs = { 'Wavelength'   : np.linspace(400e-9, 1000e-9, 500),
                     'Polarization' : [0]}

    scatSet    = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

    sourceSet  = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = None)

    Eff = Experiment.Efficiencies(['Qsca'], AsType='pymiesim')

    Eff.Plot(x='Wavelength')


if __name__ == '__main__':
    run()
