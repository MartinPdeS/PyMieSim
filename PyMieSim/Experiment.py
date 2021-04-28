#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pprint           as pp
import numpy            as np
import pandas           as pd
from beartype           import beartype
from typing             import Union
from multiprocessing    import Process
from scipy.optimize     import minimize

from PyMieSim.Source    import PlaneWave
from PyMieSim.NdArray   import PMSArray, Opt5DArray
from PyMieSim.Detector  import LPmode, Photodiode
from PyMieSim.Scatterer import Sphere, WMSample
from PyMieSim.utils     import IO, ToList

from PyMieSim.DataFrame import ( ExperimentalDataFrame,
                                 S1S2DataFrame,
                                 EfficiencesDF,
                                 ExperimentDF)

from PyMieSim.Config    import ( MetricList,
                                 DetectorParamList,
                                 SourceParamList,
                                 DefaultConfig,
                                 DefaultConfigEff )


OUTPUTTYPE = ['optimizer','numpy', 'pymiesim']
EFFTYPE    = ['Qsca', 'Qext', 'Qabs', 'Qback', 'Qratio', 'g', 'Qpr']
exList  = Union[int, float, list, np.ndarray, tuple]
exfloat = Union[bool, int, float]
DetecArg = Union[LPmode, Photodiode, list, tuple]

Print = pp.PrettyPrinter(indent=4,sort_dicts=False)

class ScatSet(object):

    @beartype
    def __init__(self,
                 DiameterList  :  exList     = None,
                 IndexList     :  exList     = None,
                 nMedium       :  exfloat    = 1.0,
                 ScattererType :  str        = 'Sphere',
                 MaterialList  : list        = None,
                 Configuration : dict        = {}):

        if MaterialList:
            assert IndexList is None, IO( "You should either choose a material or the RI, not both." )
            self.Material = MaterialList


        if IndexList:
            assert MaterialList is None, IO( "You should either choose a material or the RI, not both." )
            self.Material = None


        IndexList, DiameterList, MaterialList = ToList(IndexList, DiameterList, MaterialList)

        self._Diameter, self._Index = None, None

        self.Diameter, self.Index = DiameterList, IndexList

        self.nMedium = nMedium

        self.shape = np.shape(self.Diameter) + np.shape(self.Index)



    def UpdateConfiguration(self, config):
        i = config['MaxOrder']

        if self.Material:
            config['material'] = True
            config['order']['material']     = i
            config['label']['material']     = 'material'
            config['format']['material']    = '10s'
            config['dimension']['material'] = [mat.__name__ for mat in self.Material]
            i += 1

        else:
            config['material']              = False
            config['order']['ri']           = i
            config['label']['ri']           = 'Refracive index'
            config['format']['ri']          = '.2f'
            config['dimension']['ri']       = self.Index
            i += 1

        config['order']['diameter']     = i
        config['label']['diameter']     = 'Diameter [m]'
        config['format']['diameter']    = '.1e'
        config['dimension']['diameter'] = self.Diameter
        i += 1

        config['MaxOrder'] = i

        return config


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


    def Generator(self, source=None):
        if self.Material:
            for material in self.Material:
                for diameter in self.Diameter:
                    yield Sphere(Diameter  = diameter,
                                  Source    = source,
                                  Index     = material.Index,
                                  nMedium   = 1.)

        else:
            for diameter in self.Diameter:
                for index in self.Index:
                    yield  Sphere(Diameter  = diameter,
                                  Source    = source,
                                  Index     = index,
                                  nMedium   = 1.)


class SourceSet(object):
    @beartype
    def __init__(self,
                 WavelengthList   :   exList,
                 PolarizationList :   exList = [0],
                 SourceType       :   str    = 'PlaneWave'):

        PolarizationList, WavelengthList = ToList(PolarizationList, WavelengthList)

        self._Wavelength, self._Polarization = None, None

        self.Wavelength   = WavelengthList

        self.Polarization = PolarizationList

        self.SourceType   = SourceType

        self.shape        = np.shape(self.Wavelength) + np.shape(self.Polarization)

        self.Material     = None


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

    def Generator(self, MatGen=None):
        for wavelength in self.Wavelength:

            if self.Material:
                for material in self.Material: material.counter += 1

            for polarization in self.Polarization:

                yield PlaneWave(Wavelength   = wavelength,
                                Polarization = polarization,
                                E0           = 1)


    def UpdateConfiguration(self, config):
        i = config['MaxOrder']
        MaxOrder = np.max(list(config['order'].values()))

        config['order']['wavelength']       = i
        config['label']['wavelength']       = '$Wavelength \lambda$ [m]'
        config['format']['wavelength']      = '.1e'
        config['dimension']['wavelength']   = self.Wavelength
        i += 1

        config['order']['polarization']     = i
        config['label']['polarization']     = 'Polarization [Degree]'
        config['format']['polarization']    = '.1f'
        config['dimension']['polarization'] = self.Polarization
        i += 1

        config['MaxOrder'] = i

        return config



class DetectorSet(object):

    @beartype
    def __init__(self, DetectorList : exList):

        self.Detector = ToList(DetectorList)


    def UpdateConfiguration(self, config):
        i = config['MaxOrder']
        config['order']['detector']     = config['MaxOrder']
        config['label']['detector']     = 'Detector'
        config['format']['detector']    = '10s'
        config['dimension']['detector'] = [Det.Name for Det in self.Detector]
        config['MaxOrder']             += 1

        return config


    def Generator(self, MatGen=None):
        for detector in self.Detector:
                yield detector


class Setup(object):

    @beartype
    def __init__(self,
                 ScattererSet : ScatSet            = None,
                 SourceSet    : SourceSet          = None,
                 DetectorSet  : DetectorSet        = None):


        config = { 'name'      : None,
                   'unit'      : None,
                   'material'  : None,
                   'order'     : {},
                   'label'     : {},
                   'format'    : {},
                   'dimension' : {},
                   'MaxOrder'  : 0}

        self.MetricList   = MetricList

        self.DetectorSet  = DetectorSet

        for nd, dectector in enumerate(self.DetectorSet.Detector):
            dectector.Name = f"Detector {nd}"

        self.SourceSet    = SourceSet

        self.ScattererSet = ScattererSet

        config = self.DetectorSet.UpdateConfiguration(config)

        config = self.SourceSet.UpdateConfiguration(config)

        config = self.ScattererSet.UpdateConfiguration(config)

        self.GetShape(config)

        self.config = config


    def Efficiencies(self, Eff='Qsca', AsType='numpy'):
        """Methode generate a Pandas Dataframe of scattering efficiencies
        (Qsca) vs. scatterer diameter vs. scatterer refractive index.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe containing Qsca vs. Wavelength, Diameter vs. Index.

        """
        self.config['name']               = 'Efficiencies'
        self.config['NameList']           = Eff
        self.config['format']['variable'] = '15s'
        self.config['label']['variable']  = 'Coupling'
        self.config['unit']               = ' [1]'
        self.config['nameList']           = Eff
        self.config['output']             = AsType
        #self.config['order']['eff']       = self.config['MaxOrder']
        self.config['shape']              = self.config['shape'] + [len(Eff)]

        Eff = ToList(Eff)

        self.AssertionType(AsType=AsType)

        Array = np.empty(self.config['size'] * len(Eff))

        if self.ScattererSet.Material: self.BindMaterial()

        i = 0
        for source in self.SourceSet.Generator():
            for scatterer in self.ScattererSet.Generator(source):
                for eff in Eff:
                    Array[i]  =  getattr(scatterer, eff)
                    i += 1;

        Array = Array.reshape(self.config['shape'], order='C')

        return self.ReturnType(Array     = Array,
                               AsType    = AsType,
                               conf      = self.config)


    def AssertionType(self, AsType=None, Eff=None):
        if AsType:
            assert AsType.lower() in OUTPUTTYPE, IO( f'Invalid type {AsType}, valid choices are {OUTPUTTYPE}' )

        if Eff:
            assert set(Eff).issubset(EFFTYPE), IO( f'Invalid efficiency {Eff}, valid choices are {EFFTYPE}' )


    def Coupling(self, AsType='numpy'):
        """Property method which return a n by m by l OptArray array, n being the
        number of detectors, m is the point evaluated for the refractive index,
        l is the nomber of point evaluted for the scatterers diameters.

        Returns
        -------
        OptArray
            Raw array of detectors coupling.

        """

        self.config['name']               = 'Coupling'
        self.config['format']['variable'] = '15s'
        self.config['label']['variable']  = 'Coupling'
        self.config['unit']               = ' Watt'
        self.config['output']             = AsType

        self.AssertionType(AsType=AsType)

        Array = np.empty(self.config['size'])

        if self.ScattererSet.Material:
            self.BindMaterial()


        i = 0
        for detector in self.DetectorSet.Generator():
            for source in self.SourceSet.Generator():
                for scatterer in self.ScattererSet.Generator(source):

                    Array[i] = detector.Coupling(Scatterer = scatterer)
                    i += 1;

        return self.ReturnType(Array     = Array.reshape(self.config['shape'], order='F'),
                               AsType    = AsType,
                               conf      = self.config)


    def BindMaterial(self):
        self.SourceSet.Material = self.ScattererSet.Material

        for mat in self.ScattererSet.Material:
            mat.Evaluate(self.SourceSet.Wavelength)


    def ReturnType(self, Array, AsType, conf):
        if AsType.lower() == 'optimizer':
            return Opt5DArray(Array)

        elif AsType.lower() == 'numpy':
            return Array

        elif AsType.lower() == 'pymiesim':
            return PMSArray(array = Array, conf = conf)


    def MakeDF(self, conf, Array):

        MI = pd.MultiIndex.from_product(list(conf['dimension'].values()),
                                        names = list(conf['dimension'].keys()))

        if conf['name'].lower() == 'efficiencies':
            return EfficiencesDF(Array.reshape([conf['size'], len(config['NameList'])]),
                                 index   = MI,
                                 columns = config['NameList'])


        elif  conf['name'].lower() == 'coupling':
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


    def Optimize(self, *args, **kwargs):
        return Optimizer(Setup = self, *args, **kwargs)


class Optimizer:

    @beartype
    def __init__(self,
                 Setup         : Setup,
                 Metric        : str,
                 Parameter     : list,
                 X0            : list,
                 WhichDetector : int,
                 MinVal        : list,
                 MaxVal        : list,
                 Optimum       : str,
                 FirstStride   : Union[float, int],
                 MaxIter       : int               = 50,
                 Tol           : Union[float, int] = 1e-10):

        assert Metric.lower() in MetricList, IO( f"Metric {Metric} not in the MetricList \n{MetricList}" )
        assert all(len(x)==len(Parameter) for x in [X0, MinVal, MaxVal ]  ), IO( f'Lenght of parameters, X0, MinVal, MaxVal not equal' )

        self.Setup           = Setup
        self.Metric          = Metric
        self.Parameters      = Parameter
        self.X0              = X0
        self.WhichDetector   = WhichDetector
        self.MinVal          = MinVal
        self.MaxVal          = MaxVal
        self.FirstStride     = FirstStride
        self.MaxIter         = MaxIter
        self.Tol             = Tol

        if Optimum.lower()   == 'maximum': self.sign = -1
        elif Optimum.lower() == 'minimum': self.sign = 1

        self.Result = self.Run()


    def ComputePenalty(self, Parameters, x, MaxVal, MinVal, factor=100):
        Penalty = 0
        for n in range(len(Parameters)):
            if MinVal[n] and x[0]< MinVal[n]:
                Penalty += np.abs( x[0]*factor );
                x[0]     = self.MinVal[n]

            if MinVal[n] and x[0]> MaxVal[n]:
                Penalty += np.abs( x[0]*factor );
                x[0]     = self.MaxVal[n]

        return Penalty


    def UpdateConfiguration(self, Parameters, x, WhichDetector):

        for n in range(len(Parameters)):
            if Parameters[n] in DetectorParamList:
                setattr(self.Setup.DetectorSet[WhichDetector], Parameters[0], x[0])

            elif Parameters[n] in SourceParamList:
                setattr(self.Setup.SourceSet.Source, Parameters[0], x[0])


    def Run(self):

        def EvalFunc(x):
            Penalty = self.ComputePenalty(self.Parameters, x, self.MaxVal, self.MinVal, factor=100)

            self.UpdateConfiguration(self.Parameters, x, self.WhichDetector)

            Array = self.Setup.Coupling(AsType='Optimizer')

            Array.DefineCostFunc(self.Metric)

            return self.sign * np.abs(Array.Cost()) + Penalty

        Minimizer = Caller(EvalFunc, ParameterName = self.Parameters)

        return minimize(fun      = Minimizer.optimize,
                        x0       = self.X0,
                        method   = 'COBYLA',
                        tol      = self.Tol,
                        options  = {'maxiter': self.MaxIter, 'rhobeg':self.FirstStride})


class Caller:
    def __init__(self, function, ParameterName: list):
        self.ParameterName = ParameterName
        self.f = function # actual objective function
        self.num_calls = 0 # how many times f has been called
        self.callback_count = 0 # number of times callback has been called, also measures iteration count
        self.list_calls_inp = [] # input of all calls
        self.list_calls_res = [] # result of all calls
        self.decreasing_list_calls_inp = [] # input of calls that resulted in decrease
        self.decreasing_list_calls_res = [] # result of calls that resulted in decrease
        self.list_callback_inp = [] # only appends inputs on callback, as such they correspond to the iterations
        self.list_callback_res = [] # only appends results on callback, as such they correspond to the iterations

    def optimize(self, x):
        """Executes the actual simulation and returns the result, while
        updating the lists too. Pass to optimizer without arguments or
        parentheses."""
        result = self.f(x) # the actual evaluation of the function
        if not self.num_calls: # first call is stored in all lists
            self.decreasing_list_calls_inp.append(x)
            self.decreasing_list_calls_res.append(result)
            self.list_callback_inp.append(x)
            self.list_callback_res.append(result)
        elif result < self.decreasing_list_calls_res[-1]:
            self.decreasing_list_calls_inp.append(x)
            self.decreasing_list_calls_res.append(result)
        self.list_calls_inp.append(x)
        self.list_calls_res.append(result)
        self.num_calls += 1


        if len(self.ParameterName) == 1:

            text = """ \
            Call Number : {0} \
            \t {1}: {2:.5e} \
            \t Cost+Penalty: {3:.10e} \
            """.format(self.num_calls,
                       self.ParameterName[0],
                       x[0],
                       result)

        if len(self.ParameterName) == 2:
            text = """ \
            Call Number : {0} \
            \t {1}: {2:.5e} \
            \t {3}: {4:.5e} \
            \t Cost+Penalty: {5:.10e} \
            """.format(self.num_calls,
                       self.ParameterName[0],
                       x[0],
                       self.ParameterName[1],
                       x[1],
                       result)

        print(text)
        return result


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
