from typing  import Union
import numpy as np
import copy  as cp

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

BaseIndex = { 'name'      : None,
              'order'     :  None,
              'label'     : 'Refractive index',
              'format'    : '.2f',
              'unit'      : ' [1]',
              'dimension' : None, }

BaseDiameter = { 'name'      : 'Diameter',
                 'order'     : None,
                 'label'     : 'Diameter',
                 'format'    : '.1e',
                 'unit'      : ' m',
                 'dimension' : None,}

nMediumDict    = cp.copy(BaseIndex); nMediumDict['name']    = 'nMedium'
IndexDict      = cp.copy(BaseIndex); IndexDict['name']      = 'Index'
CoreIndexDict  = cp.copy(BaseIndex); CoreIndexDict['name']  = 'Core Index'
ShellIndexDict = cp.copy(BaseIndex); ShellIndexDict['name'] = 'Shell Index'


DiameterDict      = cp.copy(BaseDiameter); DiameterDict['name']      = 'Diameter'
CoreDiameterDict  = cp.copy(BaseDiameter); CoreDiameterDict['name']  = 'Core diameter'
ShellWidthDict    = cp.copy(BaseDiameter); ShellWidthDict['name']    = 'Shell width'


MaterialDict = { 'name'      : 'Material',
                 'order'     :  None,
                 'label'     : 'Medium Refractive index',
                 'format'    : 's',
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


BaseConfig = { 'name'      : None,
               'unit'      : {},
               'material'  : None,
               'order'     : {},
               'label'     : {},
               'format'    : {},
               'dimension' : {},
               'MaxOrder'  : 0}


Arg2Dict = { 'Diameter'      : DiameterDict,
             'Index'         : IndexDict,
             'Material'      : MaterialDict,
             'nMedium'       : nMediumDict,
             'Polarization'  : PolarizationDict,
             'Wavelength'    : WavelengthDict,
             'CoreDiameter'  : CoreDiameterDict,
             'ShellWidth'    : ShellWidthDict,
             'CoreIndex'     : CoreIndexDict,
             'ShellIndex'    : ShellIndexDict}

OUTPUTTYPE = ['optimizer','numpy', 'pymiesim']

EFFTYPE    = ['Qsca', 'Qext', 'Qabs', 'Qback', 'Qratio', 'Qpr']

exList  = Union[int, float, list, np.ndarray, tuple]
