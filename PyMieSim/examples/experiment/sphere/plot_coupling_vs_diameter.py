"""
====================
Coupling vs diameter
====================

"""


def run():
    # %%
    # Importing the package dependencies: numpy, PyMieSim
    import numpy
    from PyMieSim.experiment import SphereSet, SourceSet, Setup, PhotodiodeSet
    from PyMieSim import measure
    from PyMieSim.materials import BK7

    # %%
    # Defining the ranging parameters for the scatterer distribution
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scatterer_set = SphereSet(
        diameter=numpy.linspace(100e-9, 3000e-9, 600),
        material=BK7,
        n_medium=1.0
    )

    # %%
    # Defining the source to be employed.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    source_set = SourceSet(
        wavelength=1200e-9,
        polarization=90,
        amplitude=1
    )

    # %%
    # Defining the detector to be employed.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    detector_set = PhotodiodeSet(
        NA=[0.15, 0.1, 0.05],
        phi_offset=-180.0,
        gamma_offset=0.0,
        sampling=600,
        polarization_filter=None
    )

    # %%
    # Defining the experiment setup
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=detector_set
    )

    # %%
    # Measuring the properties
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    data = experiment.Get(measures=measure.coupling)

    # %%
    # Plotting the results
    # ~~~~~~~~~~~~~~~~~~~~
    figure = data.plot(
        y=measure.coupling,
        x=scatterer_set.diameter,
        y_scale='linear',
        normalize=True
    )

    _ = figure.show()


if __name__ == '__main__':
    run()

# -
