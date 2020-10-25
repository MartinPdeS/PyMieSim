import numpy as np
from tqdm import tqdm
from typing import Tuple
import pandas as pd

from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling



def CouplingStat(RIList: list,
                 DiameterList: list,
                 Detector,
                 **SKwargs) -> pd.DataFrame:

    Polarization = ['Parallel', 'Perpendicular']

    MI = pd.MultiIndex.from_product([Polarization, DiameterList, RIList],
                                    names=['Polarization','Diameter','RI',])

    df = pd.DataFrame(index = MI, columns = ['Coupling'])

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress") ):

        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(diameter    = Diameter,
                               index       = RI,
                               wavelength  = SKwargs['wavelength'],
                               Meshes      = SKwargs['Meshes']
                               )


            Perp, Para = PointFieldCoupling(Detector = Detector, Source   = Source)

            df.at[('Parallel', Diameter, RI),'Coupling'] = Perp

            df.at[('Perpendicular', Diameter, RI),'Coupling'] = Para

    df.Coupling = df.Coupling.astype(float)

    df = df.assign(Mean=df.groupby(['Polarization','Diameter']).Coupling.transform('mean'))

    df = df.assign(STD=df.groupby(['Polarization','Diameter']).Coupling.transform('std'))

    df._Cmax = df.Coupling.max()

    df._Cmin = df.Coupling.min()

    df._diff = np.abs(df._Cmax - df._Cmin)

    return df


def OptimizeMono(DiameterList: list,
                 Detector):
    pass
