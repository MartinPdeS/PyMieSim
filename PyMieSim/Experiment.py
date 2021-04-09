#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd


from PyMieSim.Source import PlaneWave
from PyMieSim.Optimization import NDArray, Opt4DArray, Opt5DArray
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.DataFrame import ExperimentalDataFrame, S1S2DataFrame, EfficiencesDF, ExperimentDF
from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.LMT.Scatterer import SPHERE
from PyMieSim.Source import PlaneWave


MetricList = ["max",
              "min",
              "mean",
              "rsd+RI",
              "rsd+Diameter",
              "rsd+Polarization"
              "rsd+Wavelength"
              "rsd+Detector",
              "monotonic+RI",
              "monotonic+Diameter",
              "monotonic+Polarization",
              "monotonic+Wavelength",
              "monotonic+Detector"]


class ScatSet(object):

    def __init__(self,
                 DiameterList:    list,
                 IndexList:       list,
                 nMedium:         float   = 1.0,
                 ScattererType:   str     = 'Sphere'
                 ):

        if not isinstance(IndexList, (list, np.ndarray)): IndexList = [IndexList]

        if not isinstance(DiameterList, (list, np.ndarray)): DiameterList = [DiameterList]

        self.DiameterList, self.IndexList = DiameterList, IndexList

        self.nMedium = nMedium

        self.shape = [DiameterList, IndexList]

    @property
    def Diameter(self):
        return self.DiameterList

    @Diameter.setter
    def Diameter(self, val):
        self.DiameterList = [val]

    @property
    def Index(self):
        return self.IndexList

    @Index.setter
    def Index(self, val):
        self.IndexList = [val]

    def Generator(self, Source):
        for diameter in self.DiameterList:
            for RI in self.IndexList:
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


    @property
    def Wavelength(self):
        return self.WavelengthList

    @Wavelength.setter
    def Wavelength(self, val):
        self.WavelengthList = [val]

    @property
    def Polarization(self):
        return self.WavelengthList

    @Polarization.setter
    def Polarization(self, val):
        self.PolarizationList = [val]


    def Generator(self):
        for wavelength in self.WavelengthList:
            for polarization in self.PolarizationList:
                yield PlaneWave(Wavelength   = wavelength,
                                Polarization = polarization,
                                E0           = 1)






class Setup(object):

    def __init__(self,
                 ScattererSet: ScatSet      = None,
                 SourceSet:    SourceSet    = None,
                 DetectorSet:  tuple        = None):

        if not isinstance(DetectorSet, (list, np.ndarray)): DetectorSet = [DetectorSet]

        self.DetectorSet  = DetectorSet

        self.SourceSet    = SourceSet

        self.ScattererSet = ScattererSet

        self.DetectorSetName = []
        for nd, dectector in enumerate(self.DetectorSet):
            self.DetectorSetName.append( f"Detector {nd }" )


    def Qsca(self, AsType='numpy'):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qsca vs. Wavelength, Diameter vs. Index.

        """

        conf = {'wavelength':   self.SourceSet.WavelengthList,
                'polarization': self.SourceSet.PolarizationList,
                'diameter':     self.ScattererSet.DiameterList,
                'ri':           self.ScattererSet.IndexList}


        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                Qsca, _, _ = scat.GetEfficiencies()

                Array[i] = Qsca
                i+=1

        return self.ReturnType(Array     = Array.reshape(shape),
                               AsType    = AsType,
                               conf      = conf)


    def Qext(self, AsType='numpy'):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qext) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qext vs. Wavelength, Diameter vs. Index.

        """

        conf = {'wavelength':   self.SourceSet.WavelengthList,
                'polarization': self.SourceSet.PolarizationList,
                'diameter':     self.ScattererSet.DiameterList,
                'ri':           self.ScattererSet.IndexList}


        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                _, Qext, _ = scat.GetEfficiencies()

                Array[i] = Qext
                i+=1

        return self.ReturnType(Array     = Array.reshape(shape),
                               AsType    = AsType,
                               conf      = conf)


    def Qabs(self, AsType='numpy'):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qabs) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qabs vs. Wavelength, Diameter vs. Index.

        """

        conf = {'wavelength':   self.SourceSet.WavelengthList,
                'polarization': self.SourceSet.PolarizationList,
                'diameter':     self.ScattererSet.DiameterList,
                'ri':           self.ScattererSet.IndexList}

        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                _, _, Qabs = scat.GetEfficiencies()

                Array[i] = Qabs
                i+=1


        return self.ReturnType(Array     = Array.reshape(shape),
                               AsType    = AsType,
                               conf      = conf)


    def Coupling(self, AsType='numpy'):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """

        conf = {'detector':     self.DetectorSetName,
                'wavelength':   self.SourceSet.WavelengthList,
                'polarization': self.SourceSet.PolarizationList,
                'diameter':     self.ScattererSet.DiameterList,
                'ri':           self.ScattererSet.IndexList}

        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for nd, detector in enumerate(self.DetectorSet):
            for source in self.SourceSet.Generator():
                for scat in self.ScattererSet.Generator(Source=source):

                    Array[i] = detector.Coupling(Scatterer = scat)
                    i += 1;

        return self.ReturnType(Array     = Array.reshape(shape),
                               AsType    = AsType,
                               conf      = conf)


    def ReturnType(self, Array, AsType, conf):

        if AsType.lower() == 'optimizer':
            return Opt5DArray(Array)

        elif AsType.lower() == 'numpy':
            return Array

        elif AsType.lower() == 'ndarray':
            return NDArray(array     = Array,
                           Name      = 'Coupling',
                           paramList = list(conf.keys()))

        elif AsType.lower() == 'dataframe':
            return self.MakeDF(conf, Array, Param='Coupling')


    def MakeDF(self, conf, Array, Param):
        names = list(conf.keys())
        index = list(conf.values())

        MI = pd.MultiIndex.from_product(index, names=names)

        if Param == 'Coupling':
            return ExperimentDF(Array.flatten(), index=MI, columns=[Param])

        elif Param in ['Qsca', 'Qext', 'Qabs']:
            return EfficiencesDF(Array.flatten(), index = MI, columns = [Param])


    def GetShape(self, conf):
        shape = []
        size  = 1
        for item in conf.values():
            shape += [len(item)]
            size  *= len(item)

        return shape, size





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
