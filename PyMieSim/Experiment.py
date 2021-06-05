#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import logging
import numpy              as np
from copy                 import deepcopy
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


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

class ScatSet(Set):

    @beartype
    def __init__(self, Scatterer = None, kwargs : dict = {}):

        if all([ 'Material' in kwargs.keys(), 'index' in kwargs.keys() ] ):
            raise KeyError("You should either choose a material or the RI, not both.")

        self.kwargs      = {k: ToList(v) for k, v in kwargs.items()}

        self._Scatterer_ = Scatterer

        self._Scatterer  = Namespace(kwargs=None)

        self._Source     = None


    def Generator(self):
        Generator = GeneratorFromDict(self.kwargs)

        for kwargs in Generator:
            if self._Scatterer.kwargs == kwargs:
                yield self._Scatterer

            else:
                self._Scatterer = self._Scatterer_(**kwargs, Source = self._Source)
                self._Scatterer.kwargs = deepcopy(kwargs)
                yield self._Scatterer



class SourceSet(Set):

    @beartype
    def __init__(self, Source = None, kwargs : dict = {}):

        self.kwargs   = {k: ToList(v) for k, v in kwargs.items()}

        self._Source_ = PlaneWave

        self._Source  = Namespace(kwargs=None)


    def Generator(self, MatGen=None):
        Generator = GeneratorFromDict(self.kwargs)

        for kwargs in Generator:
            if self._Source.kwargs == kwargs:
                yield self._Source

            else:

                self._Source  = self._Source_(**kwargs)
                self._Source.kwargs = deepcopy(kwargs)
                yield self._Source



class DetectorSet(Set):

    @beartype
    def __init__(self, Detector, kwargs : dict = {}):

        self.isEmpty    = False

        self.kwargs     = {k: ToList(v) for k, v in kwargs.items()}

        self._Detector_ = Detector

        self._Detector  = Namespace(kwargs=None)


    def Generator(self):
        Generator = GeneratorFromDict(self.kwargs)

        for kwargs in Generator:

            if self._Detector.kwargs == kwargs:
                yield self._Detector

            else:
                self._Detector  = self._Detector_(**kwargs)
                self._Detector.kwargs = deepcopy(kwargs)
                yield self._Detector


class EmptyDetectorSet(set):
    def __init__(self):
        self.isEmpty = True

    def UpdateConfiguration(self, config):
        return config

    def Generator(self):
        yield 1


class Setup(object):

    @beartype
    def __init__(self,
                 ScattererSet : ScatSet                  = None,
                 SourceSet    : SourceSet                = None,
                 DetectorSet  : Union[DetectorSet, None] = EmptyDetectorSet()):

        config = deepcopy(BaseConfig)

        self.SourceSet    = SourceSet

        self.DetectorSet  = DetectorSet

        self.ScattererSet = ScattererSet

        self.SourceSet.UpdateConfiguration(config)

        self.ScattererSet.UpdateConfiguration(config)

        self.DetectorSet.UpdateConfiguration(config)

        config['order'] = {dict['name']: dict['order'] for dict in config['X'].values()}

        self.config = config


    def AssertionType(self, AsType=None, Input=None):
        if 'Coupling' in Input and self.DetectorSet.isEmpty:
            raise ValueError("No coupling can be \
            computed as no detector were employed.")

        if set(PROPTYPE).intersection(Input) and not self.DetectorSet.isEmpty:
            logging.warning('The computed scatterer properties do not depends \
            on detectors although detector have been added to the experiment.')

        if AsType:
            assert AsType in OUTPUTTYPE, f'Invalid type \
            {AsType}\, valid choices are {OUTPUTTYPE}'

        if Input:
            assert set(Input).issubset(INPUTTYPE), f'Invalid \
            efficiency {Input}, valid choices are {EFFTYPE}'


    def UpdateConfig(self, Input, AsType):

        for i, prop in enumerate(Input):
            dic                 = self.config['Y']
            dic[prop]           = deepcopy( Prop2Dict[prop.lower()] )
            dic[prop]['order']  = i

        self.GetShape(self.config)

        self.config['output'] = AsType


    def Get(self, Input='Qsca', AsType='pymiesim'):
        """Methode generate array of the givens parameters as a function of
        all independent variables.

        Returns
        -------
        :class:`PyMieSimArray`
            Dataframe containing Efficiencies vs. Wavelength, Diameter vs. Index...

        """
        Input = set( ToList(Input) )

        self.AssertionType(Input=Input, AsType=AsType)

        self.UpdateConfig(Input, AsType)

        Array = np.empty(self.config['size'])

        if 'Material' in self.ScattererSet.kwargs: self.BindMaterial()

        i = 0
        for source in self.SourceSet.Generator():
            self.ScattererSet._Source = source
            for scatterer in self.ScattererSet.Generator():
                for detector in self.DetectorSet.Generator():
                    for prop in Input:
                        if prop == 'Coupling':
                            Array[i] = detector.Coupling(scatterer)
                            i       += 1

                        else:
                            Array[i] = getattr(scatterer, prop)
                            i       += 1

        Array = Array.reshape( self.config['shape'] )

        return self.ReturnType(Array = Array, AsType = AsType)


    def BindMaterial(self):
        self.SourceSet.Material = self.ScattererSet.kwargs['Material']

        for mat in self.ScattererSet.kwargs['Material']:
            mat.Evaluate(self.SourceSet.kwargs['Wavelength'])


    def ReturnType(self, Array, AsType):

        if AsType.lower() == 'optimizer':
            return Opt5DArray(Array)

        elif AsType.lower() == 'pymiesim':
            return PMSArray(array = Array, conf = self.config)


    def GetShape(self, config):

        shape = []; size  = 1
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
