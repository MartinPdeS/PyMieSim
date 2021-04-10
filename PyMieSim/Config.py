
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
                 'order'        : {
                        'detector'     : 0,
                        'wavelength'   : 1,
                        'polarization' : 2,
                        'diameter'     : 3,
                        'ri'           : 4},

                 'label'        : {
                        'variable'     : 'Coupling [Watt]',
                        'detector'     : 'Detector',
                        'wavelength'   : '$\lambda$ [m]',
                        'polarization' : 'Polarization [Degree]',
                        'diameter'     : 'Diameter [m]',
                        'ri'           : 'Refracive index'},

       }

DefaultConfigEff = {'name'         : 'Coupling',
                    'order'        : {

                            'wavelength'   : 1,
                            'polarization' : 2,
                            'diameter'     : 3,
                            'ri'           : 4},

                    'label'        : {
                            'variable'     : 'Efficiencies',
                            'wavelength'   : '$\lambda$ [m]',
                            'polarization' : 'Polarization [Degree]',
                            'diameter'     : 'Diameter [m]',
                            'ri'           : 'Refracive index'},

       }
