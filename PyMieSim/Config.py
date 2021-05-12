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


#-------------------------------------------------------------------------------
nMediumDict =   { 'name'            : 'nmedium',
                  'label'           : 'Medium refractive index [1]',
                  'format'          : '.2f',
                  'unit'            : '[1]',
                  'order'           : None,
                  'dimension'       : None }

IndexDict =     { 'name'            : 'index',
                  'label'           : 'Refractive index [1]',
                  'format'          : '.2f',
                  'unit'            : '[1]',
                  'dimension'       :  None,
                  'order'           :  None }

CoreIndexDict = { 'name'            : 'coreindex',
                  'label'           : 'Core index [1]',
                  'format'          : '.2f',
                  'unit'            : '[1]',
                  'dimension'       :  None,
                  'order'           :  None }


ShellIndexDict = { 'name'           : 'shellindex',
                   'label'          : 'Shell index [1]',
                   'format'         : '.2f',
                   'unit'           : '[1]',
                   'dimension'      :  None,
                   'order'          :  None }

DiameterDict = { 'name'             : 'diameter',
                 'label'            : 'Diameter [m]',
                 'format'           : '.1e',
                 'unit'             : '[m]',
                 'order'            : None,
                 'dimension'        : None }

CoreDiameterDict = { 'name'         : 'corediameter',
                     'label'        : 'Core diameter [m]',
                     'format'       : '.1e',
                     'unit'         : '[m]',
                     'order'        : None,
                     'dimension'    : None }

ShellWidthDict = { 'name'           : 'shellwidth',
                   'label'          : 'Shell width [m]',
                   'format'         : '.1e',
                   'unit'           : '[m]',
                   'order'          : None,
                   'dimension'      : None }

MaterialDict = { 'name'             : 'material',
                 'label'            : 'Material refractive index [1]',
                 'format'           : '10s',
                 'unit'             : '',
                 'order'            : None,
                 'dimension'        : None }

PolarizationDict = { 'name'         : 'polarization',
                     'label'        : 'Polarization [Degree]',
                     'format'       : '04.1f',
                     'unit'         : ' Degree',
                     'dimension'    : None,
                     'order'        : None }

WavelengthDict = { 'name'           : 'wavelength',
                   'label'          : r'Wavelength $\lambda$ [m]',
                   'format'         : '.1e',
                   'unit'           : 'm',
                   'order'          : None,
                   'dimension'      : None }

NADict =   { 'name'                 : 'na',
             'label'                : 'Numerical aperture (NA) [1]',
             'format'               : '.3f',
             'unit'                 : '[1]',
             'order'                : None,
             'dimension'            : None }

SamplingDict =   { 'name'           : 'sampling',
                   'label'          : 'Mesh sampling [1]',
                   'format'         : 'd',
                   'unit'           : '[1]',
                   'order'          : None,
                   'dimension'      : None }

GammaOffDict =   { 'name'           : 'gammaoffset',
                   'label'          : 'Gamma offset [degree]',
                   'format'         : '.3f',
                   'unit'           : '[degree]',
                   'order'          : None,
                   'dimension'      : None }

CouplingModeDict = { 'name'         : 'couplingmode',
                     'label'        : 'Coupling Mode',
                     'format'       : '10s',
                     'unit'         : '',
                     'order'        : None,
                     'dimension'    : None }

PhiOffDict =   { 'name'             : 'phioffset',
                 'label'            : 'Phi offset [degree]',
                 'format'           : '04.1f',
                 'unit'             : '[degree]',
                 'order'            : None,
                 'dimension'        : None }

#-------------------------------------------------------------------------------
MuScaDict = {      'type'           : u"\u03bc coefficient",
                   'name'           : 'musca',
                   'legend'         : u"\u03bc sca",
                   'label'          : "Scattering coefficient [m⁻¹]",
                   'format'         : '15s',
                   'unit'           : '[m⁻¹]' }

MuExtDict = {      'type'           : u"\u03bc coefficient",
                   'name'           : 'muext',
                   'legend'         : u"\u03bc ext",
                   'label'          : "Extinction coefficient [m⁻¹]",
                   'format'         : '15s',
                   'unit'           : '[m⁻¹]' }

MuAbsDict = {      'type'           : u"\u03bc coefficient",
                   'name'           : 'muabs',
                   'legend'         : u"\u03bc abs",
                   'label'          : "Absorption coefficient [m⁻¹]",
                   'format'         : '15s',
                   'unit'           : '[m⁻¹]' }


#-------------------------------------------------------------------------------
QscaDict = { 'type'                 : 'Efficiency',
             'name'                 : 'qsca',
             'label'                : 'Scattering efficiency [1]',
             'legend'               : "Qsca",
             'format'               : '15s',
              'unit'                : '[1]' }

QextDict = { 'type'                 : 'Efficiency',
             'name'                 : 'qext',
             'label'                : 'Extinction efficiency [1]',
             'legend'               : "Qext",
             'format'               : '15s',
              'unit'                : '[1]' }

QabsDict = { 'type'                 : 'Efficiency',
             'name'                 : 'qabs',
             'label'                : 'Absorption efficiency [1]',
             'legend'               : "Qabs",
             'format'               : '15s',
              'unit'                : '[1]' }

QratioDict = { 'type'               : 'Efficiency',
               'name'               : 'qratio',
               'label'              : 'Ratio front-back scattering efficiency [1]',
               'legend'             : "Qratio",
               'format'             : '15s',
               'unit'               : '[1]' }

QbackDict = { 'type'                : 'Efficiency',
              'name'                : 'qback',
              'label'               : 'Back-scattering efficiency [1]',
              'legend'              : "Qback",
              'format'              : '15s',
              'unit'                : '[1]' }

QprDict = { 'type'                  : 'Efficiency',
            'name'                  : 'qpr',
            'label'                 : 'Radiation pressure efficiency [1]',
            'legend'                : "Qpr",
            'format'                : '15s',
            'unit'                  : '[1]' }


#-------------------------------------------------------------------------------
CscaDict = { 'type'                 : 'Cross-section',
             'name'                 : 'csca',
             'label'                : 'Scattering cross-section [m²]',
             'legend'               : "Qsca",
             'format'               : '15s',
              'unit'                : '[m²]' }

CextDict = { 'type'                 : 'Cross-section',
             'name'                 : 'cext',
             'label'                : 'Extinction cross-section [m²]',
             'legend'               : "Qext",
             'format'               : '15s',
              'unit'                : '[m²]' }

CabsDict = { 'type'                 : 'Cross-section',
             'name'                 : 'cabs',
             'label'                : 'Absorption cross-section [m²]',
             'legend'               : "Qsca",
             'format'               : '15s',
              'unit'                : '[m²]' }

CratioDict = { 'type'               : 'Cross-section',
               'name'               : 'cratio',
               'label'              : 'Ratio front-back scattering cross-section [m²]',
               'legend'             : "Qratio",
               'format'             : '15s',
               'unit'               : '[m²]' }

CbackDict = { 'type'                : 'Cross-section',
              'name'                : 'cback',
              'label'               : 'Back-scattering cross-section [m²]',
              'legend'              : "Qback",
              'format'              : '15s',
              'unit'                : '[m²]' }

CprDict = { 'type'                  : 'Cross-section',
            'name'                  : 'cpr',
            'label'                 : 'Radiation pressure cross-section [m²]',
            'legend'                : "Qpr",
            'format'                : '15s',
            'unit'                  : '[m²]' }


#-------------------------------------------------------------------------------
gDict = { 'type'                    : 'Anisotropy factor',
          'name'                    : 'g',
          'label'                   : 'Anisotropy factor g = <cos(theta)>',
          'legend'                  : 'g',
          'format'                  : '5s',
          'unit'                    : '' }

CouplingDict = { 'name'             : 'coupling',
                 'label'            : 'Coupling [Watt]',
                 'legend'           : 'Coupling',
                 'format'           : '15s',
                 'unit'             :  '[Watt]',
                 'type'             : 'Coupling [watt]' }


#-------------------------------------------------------------------------------
BaseConfig = { 'name'               : None,
               'unit'               : {},
               'material'           : None,
               'order'              : {},
               'label'              : {},
               'format'             : {},
               'dimension'          : {},
               'MaxOrder'           : 0,
               'X'                  : {},
               'Y'                  : {}}

Arg2Dict = { 'Diameter'             : DiameterDict,
             'Index'                : IndexDict,
             'Material'             : MaterialDict,
             'nMedium'              : nMediumDict,
             'Polarization'         : PolarizationDict,
             'Wavelength'           : WavelengthDict,
             'CoreDiameter'         : CoreDiameterDict,
             'ShellWidth'           : ShellWidthDict,
             'CoreIndex'            : CoreIndexDict,
             'ShellIndex'           : ShellIndexDict,
             'NA'                   : NADict,
             'PhiOffset'            : PhiOffDict,
             'GammaOffset'          : GammaOffDict,
             'Sampling'             : SamplingDict,
             'CouplingMode'         : CouplingModeDict}

Prop2Dict = { 'musca'               : MuScaDict,
              'muext'               : MuExtDict,
              'muabs'               : MuAbsDict,
              'qsca'                : QscaDict,
              'qext'                : QextDict,
              'qabs'                : QabsDict,
              'qback'               : QbackDict,
              'qratio'              : QratioDict,
              'qpr'                 : QprDict,
              'csca'                : CscaDict,
              'cext'                : CextDict,
              'cabs'                : CabsDict,
              'cback'               : CbackDict,
              'cratio'              : CratioDict,
              'cpr'                 : CprDict,
              'g'                   : gDict,
              'coupling'            : CouplingDict  }

OUTPUTTYPE = set( ['optimizer', 'pymiesim'] )

EFFTYPE    = set( ['Qsca', 'Qext', 'Qabs', 'Qback', 'Qratio', 'Qpr'] )

CROSSTYPE  = set( ['Csca', 'Cext', 'Cabs', 'Cback', 'Cratio', 'Cpr'] )

OTHERTYPE  = set( ['g'] )

MUTYPE     = set( ['MuSca', 'MuExt', 'MuAbs'] )

PROPTYPE   = set().union(EFFTYPE, CROSSTYPE, OTHERTYPE, MUTYPE )

INPUTTYPE  = PROPTYPE.union( set( ['Coupling'] ) )

exList  = Union[int, float, list, np.ndarray, tuple]
