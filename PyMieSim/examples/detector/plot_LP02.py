"""

LP02 Mode detector
==================

"""


def run():
    from PyMieSim.detector import LPmode

    detector = LPmode(
        mode_number="LP02",
        rotation=0.,
        sampling=500,
        NA=0.6,
        gamma_offset=0,
        phi_offset=40,
        coupling_mode='Point'
    )

    figure = detector.plot()

    figure.show()


if __name__ == '__main__':
    run()

# -
