import numpy as np
from tqdm import tqdm

from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling

def OptimizeRI(RIList: list,
               DiameterList: list,
               Detector,
               **SKwargs):

    Coupling = np.empty( [ len(RIList), len(DiameterList) ] )


    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress") ):

        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(diameter    = Diameter,
                               wavelength  = SKwargs['wavelength'],
                               index       = RI,
                               npts        = SKwargs['npts'],
                               Meshes      = SKwargs['Meshes']
                               )

            Coupling[nr, nd] = PointFieldCoupling(Detector = Detector,
                                                  Source   = Source,
                                                  Field    = 'Parallel')

    STD = np.std(Coupling, axis=0)

    return STD
