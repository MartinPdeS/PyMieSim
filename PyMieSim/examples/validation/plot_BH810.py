"""
============================
Bohren-Huffman (figure~8.10)
============================

"""


def run():
    # %%
    # Importing the dependencies: numpy, matplotlib, PyMieSim
    import numpy
    import matplotlib.pyplot as plt

    from PyMieSim.tools.directories import validation_data_path
    from PyMieSim.source import PlaneWave
    from PyMieSim.scatterer import Cylinder

    theoretical = numpy.genfromtxt(f"{validation_data_path}/Figure810BH.csv", delimiter=',')

    x = theoretical[:, 0]
    y = theoretical[:, 1]

    source = PlaneWave(
        wavelength=470e-9,
        polarization=90,
        amplitude=1
    )

    scatterer_set = Cylinder(
        diameter=3000e-9,
        source=source,
        index=1.0 + 0.07j,
        n_medium=1.0
    )

    S1S2 = scatterer_set.get_s1s2(Phi=x)
    Data = (numpy.abs(S1S2.S1)**2 + numpy.abs(S1S2.S2)**2) * (0.5 / (numpy.pi * source.k))**(1 / 4)

    error = numpy.abs((Data - y) / y)

    plt.figure(figsize=(8, 4))
    plt.plot(S1S2.phi, Data, 'C1-', linewidth=3, label='PyMieSim')

    plt.plot(x, y, 'k--', linewidth=1, label='B&H [8.10]')

    plt.xlabel('scatterer_settering angle [degree]')
    plt.ylabel('Phase function')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    assert numpy.all(error < 1), 'Error: mismatch on BH_8.10 calculation occuring'


if __name__ == '__main__':
    run()

# -
