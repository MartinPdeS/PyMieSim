import numpy as np



unitList = {-6: "a",
            -5: "f",
            -4: "p",
            -3: "n",
            -2: u"\u03bc",
            -1: "m",
            +0: " ",
            +1: "K",
            +2: "M",
            +3: "G",
            +4: "T",
            +5: "P"}



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

            return f"{float(self)} {self.unit}  ({x:.2e} {unit}{self.unit})"


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


class m_3(float, UnitRep):
    """
    Class representing unit of Power [Watt]
    """
    def __new__(self, value):
        return float.__new__(self, value)

    def __init__(self, x):
        float.__init__(x)
        self.PowerFactor = -3
        self.unit = 'm⁻³'


class m_1(float, UnitRep):
    """
    Class representing unit of Power [Watt]
    """
    def __new__(self, value):
        return float.__new__(self, value)

    def __init__(self, x):
        float.__init__(x)
        self.PowerFactor = -1
        self.unit = 'm⁻¹'



# -
