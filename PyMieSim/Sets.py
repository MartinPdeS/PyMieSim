#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd


from PyMieSim.Source import PlaneWave
from PyMieSim.Optimization import OptArray
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.DataFrame import ExperimentalDataFrame, S1S2DataFrame, QscaDataFrame
from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.LMT.Scatterer import SPHERE
from PyMieSim.Source import PlaneWave

class ScatSet(object):

    def __init__(self,
                 DiameterList:    list,
                 RIList:          list,
                 nMedium:         float   = 1.0,
                 ScattererType:   str     = 'Sphere'
                 ):

        if not isinstance(RIList, (list, np.ndarray)): RIList = [RIList]

        if not isinstance(DiameterList, (list, np.ndarray)): DiameterList = [DiameterList]

        self.DiameterList, self.RIList = DiameterList, RIList

        self.nMedium = nMedium

        self.shape = [DiameterList, RIList]


    def Generator(self, Source):
        for diameter in self.DiameterList:
            for RI in self.RIList:
                yield Sphere(Diameter  = diameter,
                             Source    = Source,
                             Index     = RI,
                             nMedium   = self.nMedium,
                             MuSphere  = 1.0,
                             MuMedium  = 1.0)



class SourceSet(object):

    def __init__(self,
                 WavelengthList:      list,
                 PolarizationList:    list = [0],
                 SourceType:          str  = 'PlaneWave'):

        if not isinstance(WavelengthList, (list, np.ndarray)): WavelengthList = [WavelengthList]

        if not isinstance(PolarizationList, (list, np.ndarray)): PolarizationList = [PolarizationList]

        self.WavelengthList = WavelengthList

        self.PolarizationList = PolarizationList

        self.SourceType = SourceType

        self.shape = [WavelengthList, PolarizationList]


    def Generator(self):
        for wavelength in self.WavelengthList:
            for polarization in self.PolarizationList:
                yield PlaneWave(Wavelength   = wavelength,
                                Polarization = polarization,
                                E0           = 1)


class Experiment(object):

    def __init__(self,
                 ScattererSet: ScatSet      = None,
                 SourceSet:    SourceSet    = None,
                 DetectorSet:  tuple        = None):

        if not isinstance(DetectorSet, (list, np.ndarray)): DetectorSet = [DetectorSet]

        self.DetectorSet  = DetectorSet

        self.SourceSet    = SourceSet

        self.ScattererSet = ScattererSet


    def Efficiences(self):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca, Qext, Qabs) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qsca vs. Wavelength, Diameter vs. Index.

        """

        conf = [self.SourceSet.WavelengthList,
                self.SourceSet.PolarizationList,
                self.ScattererSet.DiameterList,
                self.ScattererSet.RIList]

        MI = pd.MultiIndex.from_product(conf, names=['Wavelength', 'Polarization', 'Diameter', 'RI'])

        df = pd.DataFrame(index = MI, columns = ['Qsca', 'Qext', 'Qabs'])

        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                Qsca, Qext, Qabs = scat.GetEfficiencies()

                df.loc[(source.Wavelength, source.Polarization.Degree, scat.Diameter, scat.Index)] = Qsca, Qext, Qabs

        return df


    def Qsca(self):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qsca vs. Wavelength, Diameter vs. Index.

        """

        conf = [self.SourceSet.WavelengthList,
                self.SourceSet.PolarizationList,
                self.ScattererSet.DiameterList,
                self.ScattererSet.RIList]

        MI = pd.MultiIndex.from_product(conf, names=['Wavelength', 'Polarization', 'Diameter', 'RI'])

        df = pd.DataFrame(index = MI, columns = ['Qsca'])

        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                Qsca, _, _ = scat.GetEfficiencies()

                df.loc[(source.Wavelength, source.Polarization.Degree, scat.Diameter, scat.Index)] = Qsca

        return df


    def Qext(self):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qext) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qext vs. Wavelength, Diameter vs. Index.

        """

        conf = [self.SourceSet.WavelengthList,
                self.SourceSet.PolarizationList,
                self.ScattererSet.DiameterList,
                self.ScattererSet.RIList]

        MI = pd.MultiIndex.from_product(conf, names=['Wavelength', 'Polarization', 'Diameter', 'RI'])

        df = pd.DataFrame(index = MI, columns = ['Qext'])

        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                _, Qext, _ = scat.GetEfficiencies()

                df.loc[(source.Wavelength, source.Polarization.Degree, scat.Diameter, scat.Index)] = Qext

        return df


    def Qabs(self):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qabs) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qabs vs. Wavelength, Diameter vs. Index.

        """

        conf = [self.SourceSet.WavelengthList,
                self.SourceSet.PolarizationList,
                self.ScattererSet.DiameterList,
                self.ScattererSet.RIList]

        MI = pd.MultiIndex.from_product(conf, names=['Wavelength', 'Polarization', 'Diameter', 'RI'])

        df = pd.DataFrame(index = MI, columns = ['Qabs'])

        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                _, _, Qabs = scat.GetEfficiencies()

                df.loc[(source.Wavelength,
                        source.Polarization.Degree,
                        scat.Diameter,
                        scat.Index)] = Qabs

        return df


    @property
    def Coupling(self):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """

        conf = [self.DetectorSet,
                self.SourceSet.WavelengthList,
                self.SourceSet.PolarizationList,
                self.ScattererSet.DiameterList,
                self.ScattererSet.RIList]

        MI = pd.MultiIndex.from_product(conf, names=['Detector',
                                                     'Wavelength',
                                                     'Polarization',
                                                     'Diameter',
                                                     'RI'])

        df = pd.DataFrame(index = MI, columns = ['Coupling'])

        for nd, detector in enumerate(self.DetectorSet):
            for source in self.SourceSet.Generator():
                for scat in self.ScattererSet.Generator(Source=source):

                    df.loc[(nd,
                            source.Wavelength,
                            source.Polarization.Degree,
                            scat.Diameter,
                            scat.Index)] = detector.Coupling(Scatterer = scat)

        return df




class ScattererSet(object):
    """Small class which represents the set of scatterer considered for analysis.

    Parameters
    ----------
    DiameterList : :class:`list`
        List of diameters to be considered for the single scatterers.
    RIList : :class:`list`
        List of refractive index to be considered for the single scatterers.
    Source : PlaneWave
        Source <class> representing the illumination impinging the scatterers.


    """

    def __init__(self,
                 DiameterList:    list,
                 RIList:          list,
                 Source:          PlaneWave,
                 nMedium:         float   = 1.0
                 ):

        if not isinstance(RIList, (list, np.ndarray)): RIList = [RIList]

        if not isinstance(DiameterList, (list, np.ndarray)): DiameterList = [DiameterList]

        self.DiameterList, self.RIList = DiameterList, RIList

        self.nMedium = nMedium

        self.Source = Source


    def Qsca(self):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing Qsca vs. Diameter vs. Index.

        """

        MI = pd.MultiIndex.from_product([self.DiameterList, self.RIList], names=['Diameter', 'RI'])

        df = QscaDataFrame(index = MI, columns = ['Qsca'])

        for nr, RI in enumerate(self.RIList):

            for nd, Diameter in enumerate(self.DiameterList):
                SizeParam =  2 * np.pi * Diameter/self.Source.Wavelength

                Scat = SPHERE(Index        = RI,
                              Diameter     = Diameter,
                              Wavelength   = self.Source.Wavelength,
                              nMedium      = self.nMedium,
                              Polarization = self.Source.Polarization.Radian,
                              E0           = self.Source.E0)

                Qsca, Qext, Qabs = Scat.Efficiencies

                df.loc[(Diameter, RI),'Qsca'] = np.abs(Qsca)

        return df


    def S1S2(self, num=201):
        """Methode generate a Pandas Dataframe of S1 and S2 parameter as a
        of the scattering angle, this for the differents scatterer diameter
        and refractive index.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing S1 and S2 vs. phi angle vs Diameter vs. Index.

        """
        Angle = np.linspace(0,2*np.pi,num)

        MI = pd.MultiIndex.from_product([self.DiameterList, self.RIList, Angle],
                                        names=['Diameter', 'RI', 'Angle'])

        df = S1S2DataFrame(index = MI, columns = ['S1', 'S2'])

        for nr, RI in enumerate(self.RIList):
            for nd, Diameter in enumerate(self.DiameterList):
                SizeParam =  2 * np.pi * Diameter/self.Source.Wavelength

                Scat = SPHERE(Index        = RI,
                              Diameter     = Diameter,
                              Wavelength   = self.Source.Wavelength,
                              nMedium      = self.nMedium,
                              Polarization = self.Source.Polarization.Radian,
                              E0           = self.Source.E0)

                S1, S2 = Scat.S1S2(Phi=np.linspace(0,2*np.pi,num))

                df.loc[(Diameter, RI),'S1'] = np.abs(S1S2[0])

                df.loc[(Diameter, RI),'S2'] = np.abs(S1S2[1])

        return df




class ExperimentalSet(object):
    """Small class which represents experimental setup for analysis.
    THe setup is comprised of ScattrererSet and a tuple of Detectors.

    Parameters
    ----------
    ScattererSet : <ScattererSet>
        Instance containing information about the scatterer to analyse in the experiment.
    Detectors : tuple
        Tuple containing all the detectors instance to analyses the scattering signal.

    """

    def __init__(self,
                 ScattererSet: ScattererSet = None,
                 Detectors:    tuple        = None):


        if isinstance(Detectors, (Photodiode, LPmode)): Detectors = (Detectors,)

        assert isinstance(Detectors, tuple), 'Detectors arguement must be a tuple'

        self.Detectors = Detectors

        self.ScattererSet = ScattererSet


    @property
    def Coupling(self):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """
        Array = np.empty( [len(self.Detectors), len(self.ScattererSet.RIList), len(self.ScattererSet.DiameterList) ] )

        for nDetector, Detector in enumerate(self.Detectors):

            for nIndex, RI in enumerate(self.ScattererSet.RIList):
                for nDiameter, Diameter in enumerate(self.ScattererSet.DiameterList):

                    Scat = Sphere(Diameter  = Diameter,
                                  Index     = RI,
                                  Source    = self.ScattererSet.Source)

                    Coupling = Detector.Coupling(Scatterer = Scat)

                    Array[nDetector, nIndex, nDiameter] = Coupling

        return OptArray(Array)


    def _Coupling(self, WhichDetector):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """
        Array = np.empty( [len(self.ScattererSet.RIList), len(self.ScattererSet.DiameterList) ] )

        for nIndex, RI in enumerate(self.ScattererSet.RIList):
            for nDiameter, Diameter in enumerate(self.ScattererSet.DiameterList):

                Scat = Sphere(Diameter  = Diameter,
                              Index     = RI,
                              Source    = self.ScattererSet.Source)

                Coupling = self.Detectors[WhichDetector].Coupling(Scatterer = Scat)

                Array[nIndex, nDiameter] = Coupling

        return OptArray(Array)

    @property
    def DataFrame(self):
        """Property method which return pandas.DataFrame of the scattering-
        detector coupling for the different diameter and refracive index
        evaluated.

        Returns
        -------
        :class:`pd.DataFrame`
            DataFrame of detectors coupling.

        """
        MI = pd.MultiIndex.from_product([range(len(self.Detectors)), self.ScattererSet.DiameterList, self.ScattererSet.RIList],
                                        names=['Detectors','Diameter','RI',])


        df = ExperimentalDataFrame(index = MI, columns = ['Coupling'])

        df.attrs['Detectors'] = self.Detectors

        for nr, RI in enumerate( self.ScattererSet.RIList ):

            for nd, Diameter in enumerate(self.ScattererSet.DiameterList):

                for nDetector, Detector in enumerate(self.Detectors):

                    Scat = Sphere(Diameter    = Diameter,
                                  Index       = RI,
                                   Source      = self.ScattererSet.Source)

                    Coupling = Detector.Coupling(Scatterer = Scat)

                    df.at[(nDetector, Diameter, RI),'Coupling'] = Coupling

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
                 Detector:        Photodiode,
                 Source:          PlaneWave,
                 Npts:            int = 201,
                 ):

        self.gList, self.LcList = gList, LcList

        self.D = D; self.Nc = Nc

        self.Detector, self.Source = Detector, Source



    @property
    def DataFrame(self):
        """Property method which return pandas.DataFrame of the scattering-
        detector coupling for the different diameter and refracive index
        evaluated.

        Returns
        -------
        :class:`pd.DataFrame`
            DataFrame of detectors coupling.

        """
        MI = pd.MultiIndex.from_product([range(len(self.Detectors)), self.ScattererSet.DiameterList, self.ScattererSet.RIList],
                                        names=['Detectors','Diameter','RI',])


        df = ExperimentalDataFrame(index = MI, columns = ['Coupling'])

        df.attrs['Detectors'] = self.Detectors

        for nr, RI in enumerate( self.ScattererSet.RIList ):

            for nd, Diameter in enumerate(self.ScattererSet.DiameterList):

                for nDetector, Detector in enumerate(self.Detectors):

                    Scat = Sample(g           = g,
                                  lc          = lc,
                                  D           = self.D,
                                  Nc          = self.Nc,
                                  Source      = LightSource,
                                  Meshes      = self.Detector.Meshes)

                    Coupling = Detector.Coupling(Scatterer = Scat)

                    df.at[(nDetector, Diameter, RI),'Coupling'] = Coupling

        df.Coupling = df.Coupling.astype(float)

        df['Mean'] = df.groupby(['Detectors','Diameter']).Coupling.transform('mean')

        df['STD'] = df.groupby(['Detectors','Diameter']).Coupling.transform('std')

        return df




    @property
    def Coupling(self):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """
        temp = np.empty( [len(self.Detectors), len(self.ScattererSet.RIList), len(self.ScattererSet.DiameterList) ] )

        for nDetector, Detector in enumerate(self.Detectors):

            for nIndex, RI in enumerate(self.ScattererSet.RIList):
                for nDiameter, Diameter in enumerate(self.ScattererSet.DiameterList):

                    Samp = Sample(g           = g,
                                  lc          = lc,
                                  D           = self.D,
                                  Nc          = self.Nc,
                                  Source      = self.Source,
                                  Meshes      = self.Detector.Meshes)

                    Coupling = Detector.Coupling(Scatterer = Samp)

                    temp[nDetector, nIndex, nDiameter] = Coupling

        return OptArray(temp)


# -
