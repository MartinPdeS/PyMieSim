import numpy as np
import matplotlib.pyplot as plt

from PyMieSim.units import ureg
from PyMieSim import single

# Fixed source
source = single.source.Gaussian(
    wavelength=750 * ureg.nanometer,
    polarization=30 * ureg.degree,
    optical_power=1 * ureg.watt,
    NA=0.3 * ureg.AU,
)

NA_values = np.linspace(0.01, 1.3, 25) * ureg.AU

# Two cases
cases = [
    ("matched  ns=1.33  nd=1.33", 1.33, 1.33),
    ("mismatch ns=1.60  nd=1.33  (TIR limited)", 1.60, 1.33),
]

results = {}

for label, ns, nd in cases:
    scatterer = single.scatterer.Sphere(
        diameter=5 * ureg.nanometer,
        source=source,
        property=1.8 * ureg.RIU,
        medium_property=ns * ureg.RIU,
    )

    scatterer.get_stokes(distance=2 * ureg.meter, sampling=80)

    couplings = []

    for NA in NA_values:
        detector = single.detector.Photodiode(
            NA=NA,
            gamma_offset=0 * ureg.AU,
            phi_offset=0 * ureg.degree,
            sampling=5000,
            medium_refractive_index=nd * ureg.RIU,
        )
        # detector = single.detector.CoherentMode(
        #     mode_number="LP11:90",
        #     NA=NA,
        #     gamma_offset=0 * ureg.AU,
        #     phi_offset=0 * ureg.degree,
        #     sampling=5000,
        #     rotation=0 * ureg.degree,
        #     medium_refractive_index=nd * ureg.RIU,
        # )
        # detector = single.detector.IntegratingSphere(
        #     sampling=5000,
        #     medium_refractive_index=nd * ureg.RIU,
        # )
        print(detector.max_angle)

        coupling_w = detector.get_coupling(scatterer=scatterer).to("watt").magnitude
        couplings.append(float(coupling_w))

    results[label] = np.array(couplings)

# Normalized curves (shape comparison)
label_matched = cases[0][0]
label_mismatch = cases[1][0]

c_matched = results[label_matched]
c_mismatch = results[label_mismatch]

c_matched_norm = c_matched / c_matched.max()
c_mismatch_norm = c_mismatch / c_mismatch.max()

ratio = c_mismatch / c_matched

# Plot
plt.figure(figsize=(10, 8))

ax1 = plt.subplot(2, 1, 1)
ax1.plot(NA_values, c_matched_norm, marker="o", label="matched (normalized)")
ax1.plot(NA_values, c_mismatch_norm, marker="o", label="ns>nd mismatch (normalized)")
ax1.set_xlabel("Photodiode NA")
ax1.set_ylabel("Normalized coupling")
ax1.set_title("Signature plot: ns>nd case should saturate earlier than matched")
ax1.grid(True)
ax1.legend()

ax2 = plt.subplot(2, 1, 2)
ax2.plot(NA_values, ratio, marker="o")
ax2.set_xlabel("Photodiode NA")
ax2.set_ylabel("Coupling ratio (ns>nd) / (matched)")
ax2.set_title("If mismatch is modeled, this ratio drops below 1 and flattens")
ax2.grid(True)

plt.tight_layout()
plt.show()

print("Matched coupling (W):", c_matched.tolist())
print("Mismatch coupling (W):", c_mismatch.tolist())
print("Ratio:", ratio.tolist())
