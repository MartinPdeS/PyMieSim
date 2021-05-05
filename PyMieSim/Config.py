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
                 'unit'      : ' [m]',
                 'dimension' : None,}

nMediumDict    = cp.copy(BaseIndex); nMediumDict['name']    = 'nmedium'   ; nMediumDict['label']    = 'Medium refractive index'
IndexDict      = cp.copy(BaseIndex); IndexDict['name']      = 'index'     ; IndexDict['name']       = 'Refractive index'
CoreIndexDict  = cp.copy(BaseIndex); CoreIndexDict['name']  = 'coreindex' ; CoreIndexDict['label']  = 'Core index'
ShellIndexDict = cp.copy(BaseIndex); ShellIndexDict['name'] = 'shellindex'; ShellIndexDict['label'] = 'Shell index'


DiameterDict      = cp.copy(BaseDiameter); DiameterDict['name']      = 'diameter';     DiameterDict['label']      = 'Diameter'
CoreDiameterDict  = cp.copy(BaseDiameter); CoreDiameterDict['name']  = 'corediameter'; CoreDiameterDict['label']  = 'Core diameter'
ShellWidthDict    = cp.copy(BaseDiameter); ShellWidthDict['name']    = 'shellwidth';   ShellWidthDict['label']    = 'Shell width'


MaterialDict = { 'name'      : 'material',
                 'order'     :  None,
                 'label'     : 'Medium Refractive index',
                 'format'    : '10s',
                 'unit'      : '',
                 'dimension' : None, }

PolarizationDict = { 'name'      : 'polarization',
                     'order'     :  None,
                     'label'     : 'Polarization',
                     'format'    : '04.1f',
                     'unit'      : ' Degree',
                     'dimension' : None,}

WavelengthDict = { 'name'      : 'wavelength',
                   'order'     :  None,
                   'label'     : r'Wavelength $\lambda$',
                   'format'    : '.1e',
                   'unit'      : ' m',
                   'dimension' : None,}

DetectorDict = { 'name'      : 'detector',
                 'order'     :  None,
                 'label'     : 'Detector',
                 'format'    : '11s',
                 'unit'      : ' ',
                 'dimension' : None, }


EfficienciesDict = { 'name'     : 'efficiencies',
                     'namelist' : None,
                     'label'    : 'Efficiencies',
                     'format'   : '15s',
                     'unit'     : ' [1]' }


CouplingDict = { 'name'   : 'coupling',
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


Variable2Unit = { 'Coupling'  : '[Watt]',
                  'Qsca'      : '[1]',
                  'Qext'      : '[1]',
                  'Qabs'      : '[1]',
                  'Qback'     : '[1]',
                  'Qratio'    : '[1]' }

Variable2Label = { 'Coupling'  : 'Coupling',
                   'Qsca'      : 'Qsca',
                   'Qext'      : 'Qext',
                   'Qabs'      : 'Qabs',
                   'Qback'     : 'Qback',
                   'Qratio'    : 'Qratio' }



OUTPUTTYPE = ['optimizer','numpy', 'pymiesim']

EFFTYPE    = ['Qsca', 'Qext', 'Qabs', 'Qback', 'Qratio', 'Qpr']

exList  = Union[int, float, list, np.ndarray, tuple]
