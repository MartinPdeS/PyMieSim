import numpy as np
import pandas as pd
from typing import Tuple, Union

from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Optimizer import OptArray
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.cpp.S1S2 import GetS1S2
from PyMieCoupling.classes.DataFrame import ExperimentalDataFrame, ScattererDataFrame
from PyMieCoupling.classes.Scattering import Scatterer, WMSample



class ScattererSet(object):

    def __init__(self,
                 DiameterList:    list,
                 RIList:          list,
                 Source:          Source):

        if not isinstance(RIList, (list, np.ndarray)): RIList = [RIList]

        if not isinstance(DiameterList, (list, np.ndarray)): DiameterList = [DiameterList]

        self.DiameterList, self.RIList = DiameterList, RIList

        self.Source = Source


    @property
    def S1S2(self, num=201):

        Angle = np.linspace(0,2*np.pi,num)

        MI = pd.MultiIndex.from_product([self.DiameterList, self.RIList, Angle],
                                        names=['Diameter', 'RI', 'Angle'])

        df = ScattererDataFrame(index = MI, columns = ['S1', 'S2'])

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):
                SizeParam =  2 * np.pi * Diameter/self.Source.Wavelength

                S1S2 = GetS1S2(RI, SizeParam, np.linspace(0,2*np.pi,num));

                df.loc[(Diameter, RI),'S1'] = np.abs(S1S2[0])

                df.loc[(Diameter, RI),'S2'] = np.abs(S1S2[1])

        return df




class ExperimentalSet(object):

    def __init__(self,
                 ScattererSet: ScattererSet = None,
                 Detectors:    list         = None):

        if not isinstance(Detectors, (list, np.ndarray)): Detectors = [Detectors]

        self.ScattererSet = ScattererSet

        self.Detectors = {'Detector {}'.format(n): detector for n, detector in enumerate(Detectors) }

        self._Coupling = None


    @property
    def Coupling(self):

        temp = np.empty( [len(self.Detectors.keys()), len(self.ScattererSet.RIList), len(self.ScattererSet.DiameterList) ] )

        for nDetector, (DetectorName, Detector) in enumerate(self.Detectors.items()):

            for nIndex, RI in enumerate(self.ScattererSet.RIList):
                for nDiameter, Diameter in enumerate(self.ScattererSet.DiameterList):

                    Scat = Scatterer(Diameter  = Diameter,
                                     Index     = RI,
                                     Source    = self.ScattererSet.Source)

                    Coupling = Detector.Coupling(Scatterer = Scat)

                    temp[nDetector, nIndex, nDiameter] = Coupling

        return OptArray(temp)


    @property
    def DataFrame(self):

        MI = pd.MultiIndex.from_product([list(self.Detectors.keys()), self.ScattererSet.DiameterList, self.ScattererSet.RIList],
                                        names=['Detectors','Diameter','RI',])


        df = ExperimentalDataFrame(index = MI, columns = ['Coupling'])

        df.attrs['Detectors'] = self.Detectors

        for nr, RI in enumerate( self.ScattererSet.RIList ):

            for nd, Diameter in enumerate(self.ScattererSet.DiameterList):

                for DetectorName, Detector in self.Detectors.items():

                    Scat = Scatterer(Diameter    = Diameter,
                                     Index       = RI,
                                     Source      = self.ScattererSet.Source)

                    Coupling = Detector.Coupling(Scatterer = Scat)

                    df.at[(DetectorName, Diameter, RI),'Coupling'] = Coupling

        df.Coupling = df.Coupling.astype(float)

        df['Mean'] = df.groupby(['Detectors','Diameter']).Coupling.transform('mean')

        df['STD'] = df.groupby(['Detectors','Diameter']).Coupling.transform('std')

        return df




class SampleSet(object):

    def __init__(self,
                 gList:           list,
                 LcList:          list,
                 D:               float,
                 Nc:              float,
                 Detector:        Union[LPmode, Photodiode],
                 Source:          Source,
                 Mode:            str = 'Centered',
                 Npts:            int = 201,
                 ):

        self.gList, self.LcList = gList, LcList

        self.D = D; self.Nc = Nc

        self.Source = Source

        self.Detector, self.Source, self.Mode = Detector, Source, Mode

        self.Coupling = np.empty( [len(self.LcList), len(self.gList)] )

        maxPhi, minPhi = np.max(Detector.Meshes.Phi.Degree), np.min(Detector.Meshes.Phi.Degree)

        self.PhiVector = np.linspace(-np.pi/2, np.pi/2,201)

        self.PhiVectorDetector = np.linspace(minPhi, maxPhi, 101)


    def GetCouplingFrame(self, Filter: list = ['None'] ):

        if not isinstance(Filter, list): Filter = [Filter]

        MI = pd.MultiIndex.from_product([Filter, self.gList, self.LcList], names=['Filter','g','lc',])

        df = DataFrame(index = MI, columns = ['Coupling'])

        df.attrs['Filter'] = Filter

        for ng, g in enumerate( self.gList ):

            for nlc, lc in enumerate(self.LcList):

                Scat = Sample(g           = g,
                              lc          = lc,
                              D           = self.D,
                              Nc          = self.Nc,
                              Source      = LightSource,
                              Meshes      = self.Detector.Meshes)

                for Polar in df.attrs['Filter']:
                    self.Detector.Filter = Polar

                    Coupling = self.Detector.Coupling(Scatterer = Scat, Mode = self.Mode)

                    df.at[(Polar, g, lc),'Coupling'] = Coupling


        df.Coupling = df.Coupling.astype(float)

        df.DetectorNane = self.Detector._name

        return df




    def GetCouplingArray(self, Filter='None'):

        temp = np.empty( [ len(self.gList), len(self.LcList) ] )


        for ng, g in enumerate( self.gList ):

            for nlc, lc in enumerate(self.LcList):

                Samp = Sample(g           = g,
                              lc          = lc,
                              D           = self.D,
                              Nc          = self.Nc,
                              Source      = self.Source,
                              Meshes      = self.Detector.Meshes)

                Coupling = self.Detector.Coupling(Scatterer = Samp, Mode = self.Mode)

                temp[ng, nlc] = Coupling

        return OptArray(temp)






# -
