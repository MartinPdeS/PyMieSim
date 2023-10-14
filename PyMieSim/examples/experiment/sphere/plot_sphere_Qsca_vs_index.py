"""
==================
Mean Qsca vs index
==================

"""


def run():
    # %%
    # Importing the package dependencies: numpy, PyMieSim
    import numpy as np
    from PyMieSim.experiment import SphereSet, SourceSet, Setup
    from PyMieSim import measure

    # %%
    # Defining the ranging parameters for the scatterer distribution
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scatterer_set = SphereSet(
        diameter=800e-9,
        index=np.linspace(1.3, 1.9, 1500),
        n_medium=1
    )

    # %%
    # Defining the source to be employed.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    source_set = SourceSet(
        wavelength=[500e-9, 1000e-9, 1500e-9],
        polarization=30,
        amplitude=1
    )

    # %%
    # Defining the experiment setup
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set
    )

    # %%
    # Measuring the properties
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    data = experiment.Get(measures=measure.Qsca)

    # %%
    # Plotting the results
    # ~~~~~~~~~~~~~~~~~~~~
    figure = data.plot(
        y=measure.Qsca,
        x=scatterer_set.index
    )

    _ = figure.show()


if __name__ == '__main__':
    run()

# -
