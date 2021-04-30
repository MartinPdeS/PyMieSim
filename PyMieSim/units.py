import numpy as np



unitList = {-6: "atto-",
            -5: "femto-",
            -4: "pico-",
            -3: "nano-",
            -2: u"\u03bc",
            -1: "milli-",
            +0: " ",
            +1: "kilo-",
            +2: "mega-",
            +3: "giga-",
            +4: "tera-",
            +5: "peta-"}



class UnitRep:
    """
    Base class for unit representation.
    """
    def __str__(self):

        if self == 0:
            return f"0 m²"

        else:
            exp = np.log10(self)//(3*self.PowerFactor)

            if   exp > +5:    unit = "P"; exp = +5
            elif exp < -5:    unit = "f"; exp = -5
            else:             unit = unitList[exp];

            x = self * 10**(-3*exp*self.PowerFactor)

            return f"{x:.2e} {unit}{self.unit}"


class Area(float, UnitRep):
    """
    Class representing unit of Area [m²]
    """
    def __new__(self, value):
        return float.__new__(self, value)

    def __init__(self, x):
        float.__init__(x)
        self.PowerFactor = 2
        self.unit = 'm²'


class Power(float, UnitRep):
    """
    Class representing unit of Power [Watt]
    """
    def __new__(self, value):
        return float.__new__(self, value)

    def __init__(self, x):
        float.__init__(x)
        self.PowerFactor = 1
        self.unit = 'Watt'






# -
