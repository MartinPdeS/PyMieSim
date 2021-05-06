def run():
    import numpy as np
    from PyMieSim            import Material
    from PyMieSim.Scatterer  import ShellSphere
    from PyMieSim.Source     import PlaneWave
    from PyMieSim.Detector   import Photodiode
    from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

    Detector0 = Photodiode(NA                = 0.1,
                           Sampling          = 300,
                           GammaOffset       = 20,
                           PhiOffset         = 30,
                           CouplingMode      = 'Centered')

    scatKwargs   = { 'CoreDiameter'     : np.geomspace(10e-09, 600e-9, 500),
                     'ShellWidth'       : [200e-9, 400e-9],
                     'CoreIndex'        : [1],
                     'ShellIndex'       : [1.3],
                     'nMedium'          : 1 }

    sourceKwargs = { 'Wavelength'   : [200e-9],
                     'Polarization' : [0]}

    scatSet   = ScatSet(Scatterer = ShellSphere,  kwargs = scatKwargs )

    sourceSet = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

    DetecSet  = DetectorSet([Detector0])

    Experiment = Setup(ScattererSet = scatSet,
                       SourceSet    = sourceSet,
                       DetectorSet  = DetecSet)

    Qsca = Experiment.Get(Properties=['Qsca', 'Qback'])

    print(Qsca)

    Qsca.Plot(x='Core diameter', Scale='lin', Groupby='name')



if __name__ == '__main__':
    run()
