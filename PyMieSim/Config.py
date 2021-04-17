
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

DefaultConfig = {'name'         : 'Coupling',
                 'unit'         : 'Watt',
                 'material'     : None,
                 'order'        : {
                                    'detector'     : 0,
                                    'wavelength'   : 1,
                                    'polarization' : 2,
                                    'diameter'     : 3,
                                    'ri'           : 4,
                                    'material'     : 4},
                 'label'        : {
                                    'variable'     : 'Coupling',
                                    'detector'     : 'Detector',
                                    'wavelength'   : '$\lambda$ [m]',
                                    'polarization' : 'Polarization [Degree]',
                                    'diameter'     : 'Diameter [m]',
                                    'ri'           : 'Refracive index',
                                    'material'     : 'material'},
                 'format'        : {
                                    'variable'     : '15s',
                                    'detector'     : '10s',
                                    'wavelength'   : '.1e',
                                    'polarization' : '.1f',
                                    'diameter'     : '.1e',
                                    'ri'           : '.2f',
                                    'material'     : '10s'},
       }


DefaultConfigEff = {'name'         : 'Efficiencies',
                    'unit'         : '[1]',
                    'material'     : None,
                    'order'        : {
                                        'wavelength'   : 1,
                                        'polarization' : 2,
                                        'diameter'     : 3,
                                        'ri'           : 4,
                                        'material'     : 4},
                    'label'        : {
                                        'variable'     : 'Efficiencies',
                                        'wavelength'   : 'Wavelength $\lambda$ [m]',
                                        'polarization' : 'Polarization [Degree]',
                                        'diameter'     : 'Diameter [m]',
                                        'ri'           : 'Refracive index',
                                        'material'     : 'material'},
                     'format'        : {
                                        'variable'     : '10s',
                                        'detector'     : '10s',
                                        'wavelength'   : '.1e',
                                        'polarization' : '.1f',
                                        'diameter'     : '.1f',
                                        'ri'           : '.2f',
                                        'material'     : '10s'},
       }
