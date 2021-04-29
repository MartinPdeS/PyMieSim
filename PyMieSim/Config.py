
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


IndexDict = { 'name'     : 'Index',
              'order'    :  None,
             'label'     : 'Refractive index',
             'format'    : '.2f',
             'dimension' : None, }

DiameterDict = {'name'      : 'Diameter',
                'order'     : None,
                'label'     : 'Diameter',
                'format'    : '.1e',
                'dimension' : None,}

nMediumDict = { 'name'      : 'nMedium',
                'order'     :  None,
                'label'     : 'Medium Refractive index',
                'format'    : '.2f',
                'dimension' : None,}

MaterialDict = { 'name'      : 'Material',
                 'order'     :  None,
                 'label'     : 'Medium Refractive index',
                 'format'    : '10s',
                 'dimension' : None, }

PolarizationDict = { 'name'      : 'Polarization',
                     'order'     :  None,
                     'label'     : 'Polarization',
                     'format'    : '.1f',
                     'dimension' : None,}

WavelengthDict = { 'name'      : 'Wavelength',
                   'order'     :  None,
                   'label'     : r'Wavelength $\lambda$',
                   'format'    : '.1e',
                   'dimension' : None,}

DetectorDict = { 'name'     : 'Detector',
                 'order'    :  None,
                 'label'     : 'Detector',
                 'format'    : '15s',
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
