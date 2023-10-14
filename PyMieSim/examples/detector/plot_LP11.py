"""

LP11 Mode detector
==================

"""


def run():
    from PyMieSim.detector import LPmode

    detector = LPmode(
        mode_number="LP11",
        rotation=0.,
        sampling=300,
        NA=0.3,
        gamma_offset=0,
        phi_offset=40,
        coupling_mode='Point'
    )

    figure = detector.plot()

    figure.show()


if __name__ == '__main__':
    run()

# -
