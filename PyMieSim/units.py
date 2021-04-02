import numpy as np



unitList = {-5: "f",
            -4: "p",
            -3: "n",
            -2: u"\u03bc",
            -1: "m",
            +0: " ",
            +1: "k",
            +2: "M",
            +3: "G",
            +4: "T",
            +5: "P"}


class Area(float):

    def __repr__(self):

        factor = 2

        exp = np.log10(self)//(3*factor)

        if   exp > +5:    unit = "P"; exp = +5
        elif exp < -5:    unit = "f"; exp = -5
        else:             unit = unitList[exp];

        x = self * 10**(-3*exp*factor)

        return f"{x:.2e} {unit}m²"




class Power(float):
    """
    P = :math:`\\int_{A} I dA`
    I = :math:`\\frac{c n \\epsilon_0}{2} |E|^2`
    With:
         I : Energy density
         n  : Refractive index of the medium
         :math:`\\epsilon_0` : Vaccum permitivity
         E  : Electric field
    """

    def __repr__(self):

        factor = 1

        exp = np.log10(self)//(3*factor)

        if   exp > +5:    unit = "P"; exp = +5
        elif exp < -5:    unit = "f"; exp = -5
        else:             unit = unitList[exp];

        x = self * 10**(-3*exp*factor)

        return f"{x:.2e} {unit}Watt"






# -
