"""
======================
coupling vs wavelength
======================

"""


def run():
    # %%
    # Importing the package dependencies: numpy, PyMieSim
    import numpy as np
    from PyMieSim.experiment import CylinderSet, SourceSet, Setup, LPModeSet
    from PyMieSim import measure
    from PyMieSim.materials import BK7

    # %%
    # Defining the ranging parameters for the scatterer distribution
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Here we look at cylinder scatterers a set diameter, refractive index and medium.
    scatterer_set = CylinderSet(
        diameter=np.linspace(100e-9, 8000e-9, 5),
        material=BK7,
        n_medium=1
    )

    # %%
    # Defining the source to be employed.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    source_set = SourceSet(
        wavelength=np.linspace(950e-9, 1050e-9, 300),
        polarization=0,
        amplitude=1
    )

    # %%
    # Defining the detector to be employed.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    detector_set = LPModeSet(
        mode_number="LP11",
        NA=[0.05, 0.01],
        phi_offset=-180,
        gamma_offset=0,
        polarization_filter=None,
        sampling=300
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
        x=source_set.wavelength,
        std=scatterer_set.diameter
    )

    _ = figure.show()


if __name__ == '__main__':
    run()

# -
