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

BaseIndex = {
              'name'      : None,
              'order'     : None,
              'label'     : None,
              'format'    : '.2f',
              'unit'      : ' [1]',
              'dimension' : None,
              }

BaseDiameter = {
                 'name'      : None,
                 'order'     : None,
                 'label'     : None,
                 'format'    : '.1e',
                 'unit'      : '[m]',
                 'dimension' : None,
                 }

nMediumDict =   {
                 'name'      : 'nmedium',
                 'label'     : 'Medium refractive index [1]',
                 'format'    : '.2f',
                 'unit'      : '[1]',
                 'order'     : None,
                 'dimension' : None,
                 }


IndexDict =     {
                 'name'      : 'index',
                 'label'     : 'Refractive index [1]',
                 'format'    : '.2f',
                 'unit'      : '[1]',
                 'dimension' :  None,
                 'order'     :  None
                 }


CoreIndexDict = {
                 'name'      : 'coreindex',
                 'label'     : 'Core index [1]',
                 'format'    : '.2f',
                 'unit'      : '[1]',
                 'dimension' :  None,
                 'order'     :  None
                 }


ShellIndexDict = {
                  'name'      : 'shellindex',
                  'label'     : 'Shell index [1]',
                  'format'    : '.2f',
                  'unit'      : '[1]',
                  'dimension' :  None,
                  'order'     :  None
                  }

DiameterDict = {
                 'name'      : 'diameter',
                 'label'     : 'Diameter [m]',
                 'format'    : '.1e',
                 'unit'      : '[m]',
                 'order'     : None,
                 'dimension' : None,
                 }

CoreDiameterDict = {
                    'name'      : 'corediameter',
                    'label'     : 'Core diameter [m]',
                    'format'    : '.1e',
                    'unit'      : '[m]',
                    'order'     : None,
                    'dimension' : None,
                   }

ShellWidthDict = {
                    'name'      : 'shellwidth',
                    'label'     : 'Shell width [m]',
                    'format'    : '.1e',
                    'unit'      : '[m]',
                    'order'     : None,
                    'dimension' : None,
                   }

MaterialDict = {
                 'name'          : 'material',
                 'label'         : 'Medium Refractive index',
                 'format'        : '10s',
                 'unit'          : '',
                 'order'         : None,
                 'dimension'     : None,
                 }

PolarizationDict = {
                     'name'      : 'polarization',
                     'label'     : 'Polarization [Degree]',
                     'format'    : '04.1f',
                     'unit'      : ' Degree',
                     'dimension' : None,
                     'order'     : None,
                     }

WavelengthDict = {
                   'name'        : 'wavelength',
                   'label'       : r'Wavelength $\lambda$ [m]',
                   'format'      : '.1e',
                   'unit'        : 'm',
                   'order'       : None,
                   'dimension'   : None
                   }

DetectorDict = {
                 'name'          : 'detector',
                 'label'         : 'Detector',
                 'format'        : '11s',
                 'unit'          : '',
                 'dimension'     : None,
                 'order'         : None,
                  }

EfficienciesDict = { 'type'      : 'efficiency',
                     'name'      : 'efficiencies',
                     'label'     : 'Efficiencies [1]',
                     'format'    : '15s',
                     'unit'      : '[1]' }

OtherDict = { 'type'             : 'other',
              'name'             : 'other',
              'label'            : ' [1]',
              'format'           : '15s',
              'unit'             : '[1]' }

CrossSectionsDict = { 'type'     : 'crossection',
                      'name'     : 'crossections',
                      'label'    : 'Cross-section [m²]',
                      'format'   : '15s',
                      'unit'     : '[m²]' }

CouplingDict = { 'name'          : 'coupling',
                 'label'         : 'Coupling [Watt]',
                 'format'        : '15s',
                 'unit'          :  '[Watt]',
                 'type'          : 'coupling' }

BaseConfig = { 'name'            : None,
               'unit'            : {},
               'material'        : None,
               'order'           : {},
               'label'           : {},
               'format'          : {},
               'dimension'       : {},
               'MaxOrder'        : 0,
               'X'               : {},
               'Y'               : {}}

Arg2Dict = { 'Diameter'          : DiameterDict,
             'Index'             : IndexDict,
             'Material'          : MaterialDict,
             'nMedium'           : nMediumDict,
             'Polarization'      : PolarizationDict,
             'Wavelength'        : WavelengthDict,
             'CoreDiameter'      : CoreDiameterDict,
             'ShellWidth'        : ShellWidthDict,
             'CoreIndex'         : CoreIndexDict,
             'ShellIndex'        : ShellIndexDict}

OUTPUTTYPE = ['optimizer', 'pymiesim']

EFFTYPE    = ['Qsca', 'Qext', 'Qabs', 'Qback', 'Qratio', 'Qpr']

CROSSTYPE  = ['Csca', 'Cext', 'Cabs', 'Cback', 'Cratio', 'Cpr']

OTHERTYPE  = ['g']

exList  = Union[int, float, list, np.ndarray, tuple]
