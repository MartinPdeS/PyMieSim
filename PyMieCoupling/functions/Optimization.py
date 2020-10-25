import numpy as np
from tqdm import tqdm
from typing import Tuple
import pandas as pd

from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling



def OptimizeRI(RIList: list,
               DiameterList: list,
               Detector,
               **SKwargs) -> pd.DataFrame:

    Polarization = ['Parallel', 'Perpendicular']

    MI = pd.MultiIndex.from_product([Polarization, DiameterList, RIList],
                                    names=['Polarization','Diameter','Index',])

    df = pd.DataFrame(index = MI, columns = ['Values'])

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress") ):

        for nd, Diameter in enumerate(DiameterList):

            Source = Scatterer(diameter    = Diameter,
                               index       = RI,
                               wavelength  = SKwargs['wavelength'],
                               Meshes      = SKwargs['Meshes']
                               )


            Perp, Para = PointFieldCoupling(Detector = Detector, Source   = Source)

            df.at[('Parallel', Diameter, RI),'Values'] = Perp

            df.at[('Perpendicular', Diameter, RI),'Values'] = Para

    df.Values = df.Values.astype(float)

    df = df.assign(STD=df.groupby(['Polarization','Diameter']).Values.transform('std'))

    df.to_html('temp.html')

    return df



def OptimizeMono(DiameterList: list,
                 Detector):
    pass
