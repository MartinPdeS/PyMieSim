#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from beartype import beartype
from typing import Union
from multiprocessing import Process

from PyMieSim.Source import PlaneWave
from PyMieSim.Optimization import PMSArray, Opt5DArray
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.DataFrame import ExperimentalDataFrame, S1S2DataFrame, EfficiencesDF, ExperimentDF
from PyMieSim.Scatterer import Sphere, WMSample


OUTPUTTYPE = ['optimizer','numpy', 'ndarray', 'dataframe']
exList = Union[list, np.ndarray]
exfloat = Union[bool, int, float]
exArg = Union[float, int, list, np.ndarray]

class ScatSet(object):

    @beartype
    def __init__(self,
                 DiameterList:    exList,
                 IndexList:       exList,
                 nMedium:         exfloat    = 1.0,
                 ScattererType:   str        = 'Sphere'):

        self._Diameter, self._Index = None, None

        self.Diameter, self.Index = DiameterList, IndexList

        self.nMedium = nMedium

        self.shape = np.shape(self.Diameter) + np.shape(self.Index)


    @property
    def Diameter(self):
        return self._Diameter

    @Diameter.setter
    def Diameter(self, val):
        if not isinstance(val, (list, np.ndarray)): val = [val]
        self._Diameter = val

    @property
    def Index(self):
        return self._Index

    @Index.setter
    def Index(self, val):
        if not isinstance(val, (list, np.ndarray)): val = [val]
        self._Index = val

    def Generator(self, Source):
        for diameter in self.Diameter:
            for RI in self.Index:
                yield Sphere(Diameter  = diameter,
                             Source    = Source,
                             Index     = RI,
                             nMedium   = self.nMedium,
                             MuSphere  = 1.0,
                             MuMedium  = 1.0)



class SourceSet(object):

    @beartype
    def __init__(self,
                 WavelengthList:      exArg,
                 PolarizationList:    exArg = [0],
                 SourceType:          str   = 'PlaneWave'):


        self._Wavelength, self._Polarization = None, None

        self.Wavelength = WavelengthList

        self.Polarization = PolarizationList

        self.SourceType = SourceType

        self.shape = np.shape(self.Wavelength) + np.shape(self.Polarization)


    @property
    def Wavelength(self):
        return self._Wavelength

    @Wavelength.setter
    def Wavelength(self, val):
        if not isinstance(val, (list, np.ndarray)): val = [val]
        self._Wavelength = val

    @property
    def Polarization(self):
        return self._Polarization

    @Polarization.setter
    def Polarization(self, val):
        if not isinstance(val, (list, np.ndarray)): val = [val]
        self._Polarization = val


    def Generator(self):
        for wavelength in self.Wavelength:
            for polarization in self.Polarization:
                yield PlaneWave(Wavelength   = wavelength,
                                Polarization = polarization,
                                E0           = 1)






class Setup(object):

    @beartype
    def __init__(self,
                 ScattererSet: ScatSet            = None,
                 SourceSet:    SourceSet          = None,
                 DetectorSet:  Union[tuple, list] = None):

        if not isinstance(DetectorSet, (list, np.ndarray)): DetectorSet = [DetectorSet]

        self.DetectorSet  = DetectorSet

        self.SourceSet    = SourceSet

        self.ScattererSet = ScattererSet

        self.DetectorSetName = []
        for nd, dectector in enumerate(self.DetectorSet):
            self.DetectorSetName.append( f"Detector {nd}" )


    def Efficiencies(self, AsType='numpy'):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qsca vs. Wavelength, Diameter vs. Index.

        """

        assert AsType.lower() in OUTPUTTYPE, \
        f'Invalid type {AsType}, valid choices are {OUTPUTTYPE}'

        conf = {'name'         : 'efficiencies',
                'order'        : {
                        'wavelength'   : 0,
                        'polarization' : 1,
                        'diameter'     : 2,
                        'ri'           : 3},
                'dimension'    : {
                        'wavelength'   : self.SourceSet.Wavelength,
                        'polarization' : self.SourceSet.Polarization,
                        'diameter'     : self.ScattererSet.Diameter,
                        'ri'           : self.ScattererSet.Index}
               }

        self.GetShape(conf)

        Array = np.empty(conf['size']*3)


        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                Qsca, Qext, Qabs = scat.GetEfficiencies()
                Array[i]   = Qsca
                Array[i+1] = Qext
                Array[i+2] = Qabs
                i+=3

        return self.ReturnType(Array     = Array.reshape([3]+conf['shape']),
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


        assert AsType.lower() in OUTPUTTYPE, \
        f'Invalid type {AsType}, valid choices are {OUTPUTTYPE}'

        conf = {'name'         : 'efficiencies',
                'order'        : {
                        'detector'     : 0,
                        'wavelength'   : 1,
                        'polarization' : 2,
                        'diameter'     : 3,
                        'ri'           : 4},

                'dimension'    : {
                        'detector'     : self.DetectorSetName,
                        'wavelength'   : self.SourceSet.Wavelength,
                        'polarization' : self.SourceSet.Polarization,
                        'diameter'     : self.ScattererSet.Diameter,
                        'ri'           : self.ScattererSet.Index}
               }


        self.GetShape(conf)

        Array = np.empty(conf['size'])

        i = 0
        for nd, detector in enumerate(self.DetectorSet):
            for source in self.SourceSet.Generator():
                for scat in self.ScattererSet.Generator(Source=source):

                    Array[i] = detector.Coupling(Scatterer = scat)
                    i += 1;

        return self.ReturnType(Array     = Array.reshape(conf['shape']),
                               AsType    = AsType,
                               conf      = conf)



    def ReturnType(self, Array, AsType, conf):

        if AsType.lower() == 'optimizer':
            return Opt5DArray(Array)

        elif AsType.lower() == 'numpy':
            return Array

        elif AsType.lower() == 'ndarray':
            return PMSArray(array     = Array,
                            Name      = 'Coupling',
                            conf      = conf)

        elif AsType.lower() == 'dataframe':
            return self.MakeDF(conf, Array)


    def MakeDF(self, conf, Array):

        MI = pd.MultiIndex.from_product(list(conf['dimension'].values()),
                                        names = list(conf['dimension'].keys()))

        if conf['Name'] == 'efficiencies':
            return EfficiencesDF(Array.reshape([conf['size'],3]),
                                 index   = MI,
                                 columns = ['Qsca', 'Qext', 'Qabs'])


        elif  conf['Name'] == 'Coupling':
            return ExperimentDF(Array.flatten(),
                                index   = MI,
                                columns = ['Coupling'])



    def GetShape(self, conf):
        shape = []
        size  = 1
        for item in conf['dimension'].values():
            shape += [len(item)]
            size  *= len(item)

        conf['shape'] = shape
        conf['size']  = size


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
