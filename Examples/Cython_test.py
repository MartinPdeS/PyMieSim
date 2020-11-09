#import pyximport; pyximport.install()
import PyMieCoupling.functions
from PyMieCoupling.S1S2 import MieS1S2
import numpy as np

if __name__ == '__main__':
    MieS1S2(np.float(1.5),
            np.float(10),
            np.float(0.5))
