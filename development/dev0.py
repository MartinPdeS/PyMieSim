import numpy as np
from PyMieSim import units
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import PlaneWave

source = PlaneWave(
    wavelength=400 * units.nanometer,  # 500 nm wavelength
    amplitude=10.0 * units.volt / units.meter,  # Relative intensity unit
    polarization=0 * units.degree,  # x-polarized
)

scatterer = Sphere(
    diameter=260 * units.nanometer,
    property=1.9 * units.RIU,
    medium_property=1.0 * units.RIU,
    source=source,
)

E_field = scatterer.get_near_field(
    # resolution=1 * units.nanometer,
    field_components=["|E|"]
)


E_field.plot()
