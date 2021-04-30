matplotlib=False
mlab=True

def run():
    from PyMieSim.Scatterer import Sphere
    from PyMieSim.Source    import PlaneWave

    Source = PlaneWave(Wavelength   = 450e-9,
                      Polarization = 0,
                      E0           = 1)

    Scat = Sphere(Diameter    = 300e-9,
                 Source      = Source,
                 Index       = 1.4)

    Stokes = Scat.Stokes(Num=100)

    Stokes.Plot()


if __name__ == '__main__':
    run()
