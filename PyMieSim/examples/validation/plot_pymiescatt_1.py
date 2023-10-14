"""
=========================
PyMieSim vs PyMieScatt: 1
=========================

"""


def run():
    # %%
    # Importing the dependencies: numpy, matplotlib, PyMieSim, PyMieScatt
    import numpy
    import matplotlib.pyplot as plt

    from PyMieSim.experiment import SourceSet, SphereSet, Setup
    from PyMieSim import measure

    import PyMieScatt as ps

    diameters = numpy.geomspace(10e-9, 6e-6, 800)
    index = 1.4

    source_set = SourceSet(
        wavelength=632.8e-9,
        polarization=0,
        amplitude=1
    )

    scatterer_set = SphereSet(
        diameter=diameters,
        index=index,
        n_medium=1.
    )

    PyMieScatt_data = []
    for diameter in diameters:
        efficiencies = ps.MieQ(
            m=index,
            wavelength=source_set.wavelength.values,
            diameter=diameter,
        )

        PyMieScatt_data.append(efficiencies[1])

    PyMieScatt_data = numpy.asarray(PyMieScatt_data).squeeze()

    experiment = Setup(
        scatterer_set=scatterer_set,
        source_set=source_set,
        detector_set=None
    )

    data = experiment.Get(measures=measure.Qsca)
    PyMieSim_data = data.array.squeeze()

    plt.figure(figsize=(8, 4))
    plt.plot(
        diameters,
        PyMieSim_data,
        'C1-',
        linewidth=3,
        label='PyMieSim'
    )

    plt.plot(
        diameters,
        PyMieScatt_data,
        'k--',
        linewidth=1,
        label='PyMieScatt'
    )

    plt.xlabel(r'diameter [$\mu$m]')
    plt.ylabel('Scattering efficiency [Sphere]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    run()

# -
