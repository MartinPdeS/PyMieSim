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


class Area(float):

    def __repr__(self):
        factor = 2

        if self == 0:
            return f"0 m²"

        else:
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

        if self == 0:
            return f"0 Watt"

        else:

            exp = np.log10(self)//(3*factor)

            if   exp > +5:    unit = "P"; exp = +5
            elif exp < -5:    unit = "f"; exp = -5
            else:             unit = unitList[exp];

            x = self * 10**(-3*exp*factor)

            return f"{x:.2e} {unit}Watt"






# -
