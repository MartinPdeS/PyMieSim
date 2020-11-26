import pandas as pd
import matplotlib.pyplot as plt
from typing import Union
from tqdm import tqdm

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Scattering import Scatterer
import numpy as np




class DataFrameCPU(pd.DataFrame):

    def __init__(self,**kwargs):
        pd.DataFrame.__init__(self,**kwargs)
        self.Polarization = None


    @property
    def Parallel(self):
        return self.xs('Parallel')


    @property
    def Perpendicular(self):
        return self.xs('Perpendicular')


    def Plot(self, y, Polarization=None):


        for Polar in self.Polarization:
            self._plot(y, Polar)


    def _plot(self, y, Polarization):

        self.ax = self.xs(Polarization).unstack(1).plot(y       = y,
                                                      grid    = True,
                                                      figsize = (8,3),
                                                      title   = '[{0}: ] {1} signal'.format(self.DetectorNane, Polarization),
                                                      ylabel  = y,
                                                      xlabel  = r'Scatterer diameter [nm]')

        plt.legend(bbox_to_anchor=(1, 1), ncol=1)

        plt.subplots_adjust(right=0.8,)



def Frame(RIList:       list,
          DiameterList: list,
          Detector:     Union[LPmode, Photodiode],
          QuietMode:    bool = False,
          Polarization: list = ['NoFiltered'],
          **SKwargs):

    if not isinstance(Polarization, list):
        Polarization = [Polarization]


    MI = pd.MultiIndex.from_product([Polarization, DiameterList, RIList],
                                    names=['Polarization','Diameter','RI',])

    df = DataFrameCPU(index = MI, columns = ['Coupling'])

    df.Polarization = Polarization

    for nr, RI in enumerate( tqdm(RIList, total = len(RIList), desc ="Progress", disable = QuietMode) ):
        for nd, Diameter in enumerate(DiameterList):

            Scat = Scatterer(Diameter    = Diameter,
                             Index       = RI,
                             Source      = SKwargs['Source'],
                             Meshes      = Detector.Meshes)

            for Polar in Polarization:
                Coupling = Detector.Coupling(Scatterer = Scat, Polarization = Polar)
                df.at[(Polar, Diameter, RI),'Coupling'] = Coupling


    df.Coupling = df.Coupling.astype(float)

    df['Mean'] = df.groupby(['Polarization','Diameter']).Coupling.transform('mean')

    df['STD'] = df.groupby(['Polarization','Diameter']).Coupling.transform('std')

    df.DetectorNane = Detector._name

    return df





# -
