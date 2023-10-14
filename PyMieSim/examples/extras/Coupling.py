"""
Photodiode Coupling
===================
"""


def run():
    from PyMieSim.source import PlaneWave
    from PyMieSim.detector import LPmode
    from PyMieSim.scatterer import Sphere

    source = PlaneWave(
        wavelength=450e-9,
        polarization=0,
        amplitude=1
    )

    detector = LPmode(
        mode_number="1-1",
        sampling=600,
        NA=0.2,
        gamma_offset=180,
        phi_offset=0,
        coupling_mode='Point'
    )

    scatterer = Sphere(
        diameter=300e-9,
        source=source,
        index=1.4
    )

    coupling = detector.coupling(scatterer=scatterer)

    print(coupling)  # 6.566085549292496e-18 Watt  (6.57e-03 fWatt)


if __name__ == '__main__':
    run()

# -
