<<<<<<< HEAD:docs/source/auto_examples/Fields:Stokes.py
"""
Fields Stokes
=============
"""

def run():
=======
def run(Plot, Save):
>>>>>>> dce48fa98dd22065341ba281e909feaeec03deb5:examples/Fields:Stokes.py
    from PyMieSim.Scatterer import Sphere
    from PyMieSim.Source    import PlaneWave

    Source = PlaneWave(Wavelength   = 450e-9,
                      Polarization = 0,
                      E0           = 1)

    Scat = Sphere(Diameter    = 300e-9,
                 Source      = Source,
                 Index       = 1.4)

    Stokes = Scat.Stokes(Num=100)

    if Plot:
        print('DDDDDDDDDd')
        Stokes.Plot()

    if Save:
        from pathlib import Path
        dir = f'docs/images/{Path(__file__).stem}'
        Stokes.SaveFig(Directory=dir)

if __name__ == '__main__':
    run(Plot=True, Save=False)
