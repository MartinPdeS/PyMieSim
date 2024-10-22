from typing import Optional
from PyOptik.units import ureg
import pint as _pint
import pint_pandas as pint

_pint.set_application_registry(ureg)


# Define a list of base units to scale
BASE_UNITS = [
    'watt', 'volt', 'meter', 'second', 'liter', 'hertz', 'ohm', 'ampere'
]

# Define prefixes for scaling units
SCALES = ['nano', 'micro', 'milli', '', 'kilo', 'mega']


def initialize_registry(ureg: Optional[object] = None):
    """
    Initialize and set up a unit registry. This function also leaks
    the units into the global namespace for easy access throughout
    the module.

    Parameters
    ----------
    ureg: Optional[pint.UnitRegistry]
        A UnitRegistry object to use. If None, the default PintType.ureg will be used.
    """

    # If no unit registry is provided, use the default
    ureg = ureg or pint.PintType.ureg

    # Set up matplotlib integration and unit formatting
    ureg.setup_matplotlib()
    ureg.formatter.default_format = '~P'  # Compact format without units like 'meter'

    # Leak scaled units into the global namespace
    for unit in BASE_UNITS:
        for scale in SCALES:
            scaled_unit_name = scale + unit
            globals()[scaled_unit_name] = getattr(ureg, scaled_unit_name)

    # Leak commonly used specific units into the global namespace
    common_units = {
        'farad': ureg.farad,
        'joule': ureg.joule,
        'coulomb': ureg.coulomb,
        'power': ureg.watt.dimensionality,
        'kelvin': ureg.kelvin,
        'celsius': ureg.celsius,
        'particle': ureg.particle,
        'RIU': ureg.refractive_index_unit,
        'refractive_index_unit': ureg.refractive_index_unit,
        'degree': ureg.degree,
        'radian': ureg.radian,
        'AU': ureg.dimensionless,
        'distance': ureg.meter.dimensionality,
        'time': ureg.second.dimensionality,
        'volume': ureg.liter.dimensionality,
        'frequency': ureg.hertz.dimensionality,
        'Quantity': ureg.Quantity
    }

    # Leak the common units into the global namespace
    globals().update(common_units)

    globals()['ureg'] = ureg


initialize_registry()
