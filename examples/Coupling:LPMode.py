def run():
    from PyMieSim.Source   import PlaneWave
    from PyMieSim.Detector import LPmode


    Source = PlaneWave(Wavelength   = 450e-9,
                      Polarization = 0,
                      E0           = 0)

    Detector = LPmode(Mode         = (1, 1),
                     Rotation     = 0.,
                     Sampling     = 201,
                     NA           = 0.4,
                     GammaOffset  = 0,
                     PhiOffset    = 40,
                     CouplingMode = 'Centered')


    Detector.Plot()


if __name__ == '__main__':
    run()
