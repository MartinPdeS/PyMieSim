# Initialize a unit registry
import pint_pandas as pint

ureg = pint.PintType.ureg

ureg.setup_matplotlib()
Quantity = ureg.Quantity
ureg.default_format = '~P'


_scaled_units_str_list = [
    'watt', 'volt', 'meter', 'second', 'liter', 'hertz', 'ohm', 'ampere'
]

for _units_str in _scaled_units_str_list:
    for scale in ['nano', 'micro', 'milli', '', 'kilo', 'mega']:
        _unit = scale + _units_str
        globals()[_unit] = getattr(ureg, _unit)


joule = ureg.joule
coulomb = ureg.coulomb
power = ureg.watt.dimensionality
kelvin = ureg.kelvin
celsius = ureg.celsius
particle = ureg.particle
RIU = refractive_index_unit = ureg.refractive_index_unit
degree = ureg.degree
radian = ureg.radian
AU = ureg.dimensionless
distance = ureg.meter.dimensionality
time = ureg.second.dimensionality
volume = ureg.liter.dimensionality
frequency = ureg.hertz.dimensionality
