"""

Photodiode detector
===================

"""


def run():
    from PyMieSim.detector import Photodiode

    detector = Photodiode(
        NA=0.8,
        sampling=500,
        gamma_offset=0,
        phi_offset=0,
        polarization_filter=None
    )

    figure = detector.plot()

    figure.show()


if __name__ == '__main__':
    run()

# -
