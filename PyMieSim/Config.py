from typing  import Union
import numpy as np

MetricList = ["max",
              "min",
              "mean",
              "rsd+ri",
              "rsd+diameter",
              "rsd+polarization"
              "rsd+wavelength"
              "rsd+detector",
              "monotonic+ri",
              "monotonic+diameter",
              "monotonic+polarization",
              "monotonic+wavelength",
              "monotonic+detector"]

DetectorParamList = ['NA',
                     'PhiOffset',
                     'ThetaOffset',
                     'Filter']

SourceParamList = ['E0',
                   'Polarization',
                   'Wavelength']


IndexDict = { 'name'      : 'Index',
              'order'     :  None,
              'label'     : 'Refractive index',
              'format'    : '.2f',
              'unit'      : ' [1]',
              'dimension' : None, }

DiameterDict = {'name'      : 'Diameter',
                'order'     : None,
                'label'     : 'Diameter',
                'format'    : '.1e',
                'unit'      : ' m',
                'dimension' : None,}

nMediumDict = { 'name'      : 'nMedium',
                'order'     :  None,
                'label'     : 'Medium Refractive index',
                'format'    : '.2f',
                'unit'      : ' [1]',
                'dimension' : None,}

MaterialDict = { 'name'      : 'Material',
                 'order'     :  None,
                 'label'     : 'Medium Refractive index',
                 'format'    : '10s',
                 'unit'      : '',
                 'dimension' : None, }

PolarizationDict = { 'name'      : 'Polarization',
                     'order'     :  None,
                     'label'     : 'Polarization',
                     'format'    : '.1f',
                     'unit'      : ' Degree',
                     'dimension' : None,}

WavelengthDict = { 'name'      : 'Wavelength',
                   'order'     :  None,
                   'label'     : r'Wavelength $\lambda$',
                   'format'    : '.1e',
                   'unit'      : ' m',
                   'dimension' : None,}

DetectorDict = { 'name'     : 'Detector',
                 'order'    :  None,
                 'label'     : 'Detector',
                 'format'    : '11s',
                 'unit'      : ' ',
                 'dimension' : None, }


EfficienciesDict = { 'name'     : 'Efficiencies',
                     'namelist' : None,
                     'label'    : 'Efficiencies',
                     'format'   : '15s',
                     'unit'     : ' [1]' }


CouplingDict = { 'name'   : 'Coupling',
                 'label'  : 'Coupling',
                 'format' : '15s',
                 'unit'   :  ' Watt' }

Arg2Dict = { 'Diameter'     : DiameterDict,
             'Index'        : IndexDict,
             'Material'     : MaterialDict,
             'nMedium'      : nMediumDict,
             'Polarization' : PolarizationDict,
             'Wavelength'   : WavelengthDict}

OUTPUTTYPE = ['optimizer','numpy', 'pymiesim']

EFFTYPE    = ['Qsca', 'Qext', 'Qabs', 'Qback', 'Qratio', 'g', 'Qpr']

exList  = Union[int, float, list, np.ndarray, tuple]
