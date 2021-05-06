#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import copy
import numpy              as np
from beartype             import beartype
from multiprocessing      import Process
from scipy.optimize       import minimize

from PyMieSim.Source      import PlaneWave
from PyMieSim.NdArray     import PMSArray, Opt5DArray
from PyMieSim.Detector    import LPmode, Photodiode
from PyMieSim.Scatterer   import Sphere, WMSample
from PyMieSim.BaseClasses import Set
from PyMieSim.utils       import IO, ToList, GeneratorFromDict, MergeDict
from PyMieSim.Config      import *


class ScatSet(Set):

    @beartype
    def __init__(self, Scatterer = None, kwargs : dict = {}):

        if all([ 'Material' in kwargs.keys(), 'index' in kwargs.keys() ] ):
            raise KeyError("You should either choose a material or the RI, not both.")

        self.kwargs = {k: ToList(v) for k, v in kwargs.items()}

        self.Scatterer = Scatterer


    def Generator(self, source=None):
        Generator, order = GeneratorFromDict(self.kwargs)

        for arguments in Generator:
            for n, key in enumerate( self.kwargs.keys() ):
                order[key] = arguments[n]

            yield self.Scatterer(**order, Source = source)



class SourceSet(Set):

    @beartype
    def __init__(self, Source = None, kwargs : dict = {}):

        self.kwargs = {k: ToList(v) for k, v in kwargs.items()}


    def Generator(self, MatGen=None):
        Generator, order = GeneratorFromDict(self.kwargs)

        for arguments in Generator:
            for n, key in enumerate( self.kwargs.keys() ):
                order[key] = arguments[n]

            yield PlaneWave(**order)



class DetectorSet(Set):

    @beartype
    def __init__(self, DetectorList : exList):

        self.Detector = ToList(DetectorList)

        for nd, dectector in enumerate(self.Detector):
            dectector.Name = f"Detector {nd}"


    def UpdateConfiguration(self, config):

        i = config['MaxOrder']

        Dict              = copy.deepcopy( DetectorDict )
        Dict['order']     = i
        Dict['dimension'] = [Det.Name for Det in self.Detector]
        Dict['size']      = len(Dict['dimension'])

        MergeDict(config,Dict)
        config['X'][i] = Dict

        config['MaxOrder'] = i+1


        return config


    def Generator(self, MatGen=None):
        for detector in self.Detector:
                yield detector



class Setup(object):

    @beartype
    def __init__(self,
                 ScattererSet : ScatSet                  = None,
                 SourceSet    : SourceSet                = None,
                 DetectorSet  : Union[DetectorSet, None] = None):

        config = copy.deepcopy(BaseConfig)

        self.MetricList   = MetricList

        self.DetectorSet  = DetectorSet

        self.SourceSet    = SourceSet

        self.ScattererSet = ScattererSet

        if DetectorSet: config = self.DetectorSet.UpdateConfiguration(config)

        config = self.SourceSet.UpdateConfiguration(config)

        config = self.ScattererSet.UpdateConfiguration(config)

        config['order'] = {dict['name']: dict['order'] for dict in config['X'].values()}
        
        self.config = config


    def AssertionType(self, AsType=None, Eff=None):
        if AsType:
            assert AsType.lower() in OUTPUTTYPE, IO( f'Invalid type {AsType}, valid choices are {OUTPUTTYPE}' )

        if Eff:
            assert set(Eff).issubset(EFFTYPE), IO( f'Invalid efficiency {Eff}, valid choices are {EFFTYPE}' )


    def UpdateConfig(self, Properties, AsType):

        for i, prop in enumerate(Properties):

            if prop in EFFTYPE:
                self.config['Y'][prop] = copy.deepcopy( EfficienciesDict )
                self.config['Y'][prop]['name'] = prop
                self.config['Y'][prop]['order'] = i

            else:
                self.config['Y'][prop] = copy.deepcopy( CouplingDict )
                self.config['Y'][prop]['name'] = prop
                self.config['Y'][prop]['order'] = i

        self.GetShape(self.config)

        self.config['output'] = AsType


    def Get(self, Properties='Qsca', AsType='pymiesim'):
        """Methode generate array of the givens parameters as a function of
        all independent variables.

        Returns
        -------
        :class:`PyMieSimArray`
            Dataframe containing Efficiencies vs. Wavelength, Diameter vs. Index...

        """
        if Properties == 'all': Properties = EFFTYPE

        self.AssertionType(AsType=AsType)

        Properties = ToList(Properties)

        self.UpdateConfig(Properties, AsType)

        Array = np.empty(self.config['size'])

        if 'Material' in self.ScattererSet.kwargs: self.BindMaterial()

        i = 0
        for source in self.SourceSet.Generator():
            for scatterer in self.ScattererSet.Generator(source):
                for eff in Properties:
                    if eff == 'Coupling':
                        for detector in self.DetectorSet.Generator():
                            Array[i] = detector.Coupling(Scatterer = scatterer)
                            i       += 1;

                    else:
                        Array[i] = getattr(scatterer, eff)
                        i       += 1

        return self.ReturnType(Array     = Array.reshape( self.config['shape'] ),
                               AsType    = AsType,
                               conf      = self.config)


    def BindMaterial(self):
        self.SourceSet.Material = self.ScattererSet.kwargs['Material']

        for mat in self.ScattererSet.kwargs['Material']:
            mat.Evaluate(self.SourceSet.kwargs['Wavelength'])


    def ReturnType(self, Array, AsType, conf):

        if AsType.lower() == 'optimizer':
            return Opt5DArray(Array)

        elif AsType.lower() == 'pymiesim':
            return PMSArray(array = Array, conf = conf)


    def GetShape(self, config):

        shape = []
        size  = 1
        for key, val in config['X'].items():
            shape += [val['size']]
            size  *= val['size']

        length = len( [val['name'] for val in config['Y'].values()] )

        config['shape'] = shape + [length]
        config['size']  = size  *  length


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
