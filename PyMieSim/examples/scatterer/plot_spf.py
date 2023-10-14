"""
===============
SPF computation
===============

"""


def run():
    # %%
    # Importing the package: PyMieSim
    from PyMieSim.scatterer import Sphere
    from PyMieSim.source import PlaneWave

    # %%
    # Defining the source
    # ~~~~~~~~~~~~~~~~~~~
    source = PlaneWave(
        wavelength=500e-9,
        polarization=0,
        amplitude=1
    )

    # %%
    # Defining the scatterer
    # ~~~~~~~~~~~~~~~~~~~~~~
    scatterer = Sphere(
        diameter=1200e-9,
        source=source,
        index=1.4,
        n_medium=1.0
    )

    # %%
    # Computing the data
    # ~~~~~~~~~~~~~~~~~~
    data = scatterer.get_spf(sampling=300)

    # %%
    # Plotting the data
    # ~~~~~~~~~~~~~~~~~
    figure = data.plot()

    _ = figure.show()


if __name__ == '__main__':
    run()

# -
