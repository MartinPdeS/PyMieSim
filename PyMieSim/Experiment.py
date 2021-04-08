#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd


from PyMieSim.Source import PlaneWave
from PyMieSim.Optimization import OptArray
from PyMieSim.Detector import LPmode, Photodiode
from PyMieSim.DataFrame import ExperimentalDataFrame, S1S2DataFrame, EfficiencesDF, ExperimentDF
from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.LMT.Scatterer import SPHERE
from PyMieSim.Source import PlaneWave


MetricList = ["max",
              "min",
              "mean",
              "std+RI",
              "std+Diameter",
              "std+Polarization"
              "std+Wavelength"
              "std+Detector",
              "monotonic+RI",
              "monotonic+Diameter",
              "monotonic+Polarization",
              "monotonic+Wavelength",
              "monotonic+Detector"]


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



class Opt5DArray(np.ndarray):
    def __new__(cls, *args, **kwargs):
        this = np.array(*args, **kwargs, copy=False)
        this = np.asarray(this).view(cls)

        this.dimensions = ['Detector',
                           'Wavelength',
                           'Polarization',
                           'Diameter',
                           'Index']
        return this

    def __array_finalize__(self, obj):
        pass


    def __init__(self, arr):
        pass


    def Cost(self, arg = 'max'):

        arg = arg.lower().split('+', 2)

        if len(arg) == 1:
            if   'max' in arg:  return np.max(self)
            elif 'min' in arg:  return np.min(self)
            elif 'mean' in arg: return np.mean(self)

        if len(arg) == 2:
            if   arg[0] == 'rsd':        func = self.rsd
            elif arg[0] == 'monotonic':  func = self.Monotonic

            if   arg[1] == 'ri':           return np.mean( func(self, axis = 4) )
            elif arg[1] == 'diameter':     return np.mean( func(self, axis = 3) )
            elif arg[1] == 'polarization': return np.mean( func(self, axis = 2) )
            elif arg[1] == 'wavelength':   return np.mean( func(self, axis = 1) )
            elif arg[1] == 'detector':     return np.mean( func(self, axis = 0) )

        raise ValueError(f"Invalid metric input. \nList of metrics: {MetricList}")


    def Monotonic(self, axis):

        Grad = np.gradient(self, axis = axis)

        STD = Grad.std( axis = axis)

        return STD[0]


    def rsd(self, array, axis):
        return np.std(array, axis)/np.mean(array, axis)


    def RIMonotonic(self):

        Grad = np.gradient(self, axis = 0)

        STD = Grad.std( axis = 0)

        return STD[0]


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


    def Qsca(self, AsDataframe=True):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qsca vs. Wavelength, Diameter vs. Index.

        """

        conf = {'Wavelength':   self.SourceSet.WavelengthList,
                'Polarization': self.SourceSet.PolarizationList,
                'Diameter':     self.ScattererSet.DiameterList,
                'RI':           self.ScattererSet.RIList}

        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                Qsca, _, _ = scat.GetEfficiencies()

                Array[i] = Qsca
                i+=1

        if AsDataframe is False: return Opt5DArray(Array.reshape(shape))

        else:
            return self.MakeDF(conf, Array, Param='Qsca')


    def Qext(self, AsDataframe=True):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qext) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qext vs. Wavelength, Diameter vs. Index.

        """

        conf = {'Wavelength':   self.SourceSet.WavelengthList,
                'Polarization': self.SourceSet.PolarizationList,
                'Diameter':     self.ScattererSet.DiameterList,
                'RI':           self.ScattererSet.RIList}

        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                _, Qext, _ = scat.GetEfficiencies()

                Array[i] = Qext
                i+=1

        if AsDataframe is False: return Opt5DArray(Array.reshape(shape))

        else:
            return self.MakeDF(conf, Array, Param='Qext')


    def Qabs(self, AsDataframe=True):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qabs) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qabs vs. Wavelength, Diameter vs. Index.

        """

        conf = {'Wavelength':   self.SourceSet.WavelengthList,
                'Polarization': self.SourceSet.PolarizationList,
                'Diameter':     self.ScattererSet.DiameterList,
                'RI':           self.ScattererSet.RIList}

        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for source in self.SourceSet.Generator():
            for scat in self.ScattererSet.Generator(Source=source):
                _, _, Qabs = scat.GetEfficiencies()

                Array[i] = Qabs
                i+=1


        if AsDataframe is False: return Opt5DArray(Array.reshape(shape))

        else:
            return self.MakeDF(conf, Array, Param='Qabs')


    def Coupling(self, AsDataframe=False):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """

        conf = [self.DetectorSetName,
                self.SourceSet.WavelengthList,
                self.SourceSet.PolarizationList,
                self.ScattererSet.DiameterList,
                self.ScattererSet.RIList]

        shape, size = self.GetShape(conf)

        Array = np.empty(size)

        i = 0
        for nd, detector in enumerate(self.DetectorSet):
            for source in self.SourceSet.Generator():
                for scat in self.ScattererSet.Generator(Source=source):

                    Array[i] = detector.Coupling(Scatterer = scat)
                    i += 1;

        if AsDataframe is False: return Opt5DArray(Array.reshape(shape))

        else:
            self.MakeDF(conf, Array, Param='Coupling')


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