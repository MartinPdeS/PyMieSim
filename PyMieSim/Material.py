import numpy as np


#https://en.wikipedia.org/wiki/Sellmeier_equation
def BK7Glass(wavelength):
    B1 = 1.03961212
    B2 = 0.231792344
    B3 = 1.01046945

    C1 = 6.00069867e-3 / 1e12 #[micro M^2]
    C2 = 2.00179144e-2 / 1e12 #[micro M^2]
    C3 = 1.03560653e+2 / 1e12 #[micro M^2]

    return Getn(A          = 1,
                B          = (B1, B2, B3),
                C          = (C1, C2, C3),
                wavelength = wavelength)


#https://en.wikipedia.org/wiki/Sellmeier_equation
def FusedSilica(wavelength):
    B1 = 0.696166300
    B2 = 0.407942600
    B3 = 0.897479400

    C1 = 4.67914826e-3 / 1e12 #[micro M^2]
    C2 = 1.35120631e-2 / 1e12#[micro M^2]
    C3 = 9.79340025e+1 / 1e12#[micro M^2]

    return Getn(A          = 1,
                B          = (B1, B2, B3),
                C          = (C1, C2, C3),
                wavelength = wavelength)


#https://en.wikipedia.org/wiki/Sellmeier_equation
def AmbiantAir(wavelength):
    B1 = 4.915889e-4
    B2 = 5.368273e-5
    B3 = -1.949567e-4

    C1 = 4.352140e-3 / 1e12 #[micro M^2]
    C2 = 1.747001e-2 / 1e12 #[micro M^2]
    C3 = 4.258444e+3 / 1e12 #[micro M^2]

    return Getn(A          = 1,
                B          = (B1, B2, B3),
                C          = (C1, C2, C3),
                wavelength = wavelength)


 # https://www.researchgate.net/figure/Sellmeier-coefficients-of-As-2-S-5-and-borosilicate-glasses_tbl1_303895234
def BorosilicateGlass(wavelength):
    B1 = 0.020452
    B2 = 107.9261
    B3 = 0.0002333

    C1 = 2.1361 / 1e12 #[micro M^2]
    C2 = 0.0693 / 1e12 #[micro M^2]
    C3 = 1.7637 / 1e12 #[micro M^2]

    return Getn(A          = 1,
                B          = (B1, B2, B3),
                C          = (C1, C2, C3),
                wavelength = wavelength)



#https://www.researchgate.net/figure/Sellmeier-parameters-for-water-24_tbl1_315530939
def Water(wavelength):
    B1 = 0.5684027565
    B2 = 0.1726177391
    B3 = 0.02086189578

    C1 = 5.101829712e-3 / 1e12 #[micro M^2]
    C2 = 1.821153936e-2 / 1e12 #[micro M^2]
    C3 = 2.620722293e-2 / 1e12 #[micro M^2]

    return Getn(A          = 1,
                B          = (B1, B2, B3),
                C          = (C1, C2, C3),
                wavelength = wavelength)


def Getn(A, B, C, wavelength):
    term0 = A
    term1 = B[0]*wavelength**2/(wavelength**2 - C[0])
    term2 = B[1]*wavelength**2/(wavelength**2 - C[1])
    term3 = B[2]*wavelength**2/(wavelength**2 - C[2])
    sum = term0 + term1 + term2 + term3
    #sum = 2.6708254257402237
    print('$$$$',sum)
    print('####',sum , np.sqrt(sum))

    return np.sqrt(sum)

#-
