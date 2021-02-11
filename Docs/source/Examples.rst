Examples
========

Scatterer: S1-S2
----------------

.. code-block:: console
   :linenos:

   from PyMieCoupling.classes.Scattering import Scatterer
   from PyMieCoupling.Physics import Source

   LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

   Scat = Scatterer(Diameter    = 300e-9,
                    Source      = LightSource,
                    Index       = 1.4)


   S1S2 = Scat.S1S2(Num=100)

   S1S2.Plot()


.. image:: ../images/S1S2.png
   :width: 600

Scatterer: full far-field
-------------------------

.. code-block:: console
   :linenos:

   from PyMieCoupling.classes.Scattering import Scatterer
   from PyMieCoupling.Physics import Source

   LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

   Scat = Scatterer(Diameter    = 300e-9,
                    Source      = LightSource,
                    Index       = 1.4)


   Fields = Scat.Field(Num=100)

   Fields.Plot()

.. image:: ../images/Fields.png
   :width: 600

Scatterer: phase function
-------------------------

.. code-block:: console
   :linenos:

   from PyMieCoupling.classes.Scattering import Scatterer
   from PyMieCoupling.Physics import Source

   LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

   Scat = Scatterer(Diameter    = 400e-9,
                    Source      = LightSource,
                    Index       = 1.4)


   SPF = Scat.SPF(Num=100)

   SPF.Plot()

.. image:: ../images/SPF.png
   :width: 600

Detector: Photodiode
--------------------

.. code-block:: console
   :linenos:

   from PyMieCoupling.Physics import Source
   from PyMieCoupling.classes.Detector import Photodiode

   LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

   Detector = Photodiode(NA                = 0.8,
                         Sampling          = 1001,
                         GammaOffset       = 0,
                         PhiOffset         = 0)


   Detector.Plot()


.. image:: ../images/Photodiode.png
   :width: 600

Detector: LPMode
----------------

.. code-block:: console
   :linenos:

   from PyMieCoupling.Physics import Source
   from PyMieCoupling.classes.Detector import LPmode

   LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

   Detector = LPmode(Mode         = (1, 1,'h'),
                     Sampling     = 201,
                     NA           = 0.2,
                     GammaOffset  = 0,
                     PhiOffset    = 0,
                     CouplingMode = 'Centered')


   Detector.Plot()

.. image:: ../images/LPmode.png
   :width: 600

Coupling: Scatterer-LPMode
--------------------------

.. code-block:: console
   :linenos:

   from PyMieCoupling.Physics import Source
   from PyMieCoupling.classes.Detector import LPmode
   from PyMieCoupling.classes.Scattering import Scatterer

   LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

   Detector = LPmode(Mode         = (1, 1,'h'),
                     Sampling     = 201,
                     NA           = 0.2,
                     GammaOffset  = 0,
                     PhiOffset    = 0,
                     CouplingMode = 'Centered')


   Scat = Scatterer(Diameter    = 400e-9,
                    Source      = LightSource,
                    Index       = 1.4)

   Coupling = Detector.Coupling(Scatterer = Scat)

   print(Coupling)


Output: (2.852590820006693e-07)


ScattererSet: Qscattering
--------------------------

.. code-block:: console
   :linenos:

   import numpy as np
   from PyMieCoupling.Physics import Source
   from PyMieCoupling.classes.Sets import ScattererSet

   LightSource = Source(Wavelength   = 950e-9,
                        Polarization = 0)


   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 500e-9, 100),
                          RIList        = np.linspace(1.5, 1.5, 1).round(1),
                          Source        = LightSource)


   Qsca = ScatSet.Qsca()

   Qsca.Plot()


.. image:: ../images/Qsca.png
   :width: 600

ExperimentalSet: Coupling
----------------------------

.. code-block:: console
   :linenos:

   import numpy as np
   from PyMieCoupling.Physics import Source
   from PyMieCoupling.classes.Detector import LPmode
   from PyMieCoupling.classes.Sets import ScattererSet, ExperimentalSet

   LightSource = Source(Wavelength   = 950e-9,
                        Polarization = 0)



   Detector0 = LPmode(NA                = 0.2,
                      Sampling          = 401,
                      GammaOffset       = 0,
                      PhiOffset         = 20,
                      Mode              = (0,1),
                      CouplingMode      = 'Mean')

   Detector1 = LPmode(NA                = 0.2,
                      Sampling          = 401,
                      GammaOffset       = 0,
                      PhiOffset         = 20,
                      Mode              = (1,1),
                      CouplingMode      = 'Mean')





   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 3000e-9, 300),
                          RIList        = np.linspace(1.5, 1.5, 1).round(1),
                          Source        = LightSource)





   Set = ExperimentalSet(ScattererSet  = ScatSet,
                         Detectors     = [Detector0, Detector1])


   Data = Set.DataFrame

   Data.Plot(y='Coupling')


.. image:: ../images/ExperimentalSet.png
   :width: 600











Optimizer: NA
----------------------------

.. code-block:: console
  :linenos:





   import numpy as np
   from scipy.optimize import minimize
   from PyMieCoupling.classes.Detector import Photodiode, LPmode
   from PyMieCoupling.Physics import Source
   from PyMieCoupling.classes.Optimizer import Simulator
   from PyMieCoupling.classes.Sets import ExperimentalSet, ScattererSet

   LightSource = Source(Wavelength = 450e-9,
                        Polarization = 0,
                        Power = 1,
                        Radius = 1)



   Detector0 = Photodiode(NA                = 0.2,
                          Sampling          = 150,
                          GammaOffset       = 0,
                          PhiOffset         = 0,
                          CouplingMode      = 'Centered')

   Detector1 = LPmode(NA                = 0.2,
                      Sampling          = 150,
                      Mode              = (0,1),
                      GammaOffset       = 0,
                      PhiOffset         = 0,
                      CouplingMode      = 'Centered')


   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 3500e-9, 100),
                          RIList        = np.linspace(1.5, 1.5, 1).round(1),
                          Source        = LightSource)

   Set = ExperimentalSet(ScattererSet  = ScatSet,
                         Detectors     = (Detector0, Detector1))


   def EvalFunc(x):

       Set.Detectors[1].NA = x[0]

       return Set.Coupling.Cost('Max') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


   Minimizer = Simulator(EvalFunc, ParameterName= ['NA'])

   Result = minimize(fun      = Minimizer.simulate,
                     x0       = [0.2],
                     method   = 'COBYLA',
                     callback = Minimizer.callback,
                     tol      = 1e-5,
                     options  = {'maxiter': 10, 'rhobeg':0.1})

   print(Result)

   Set.DataFrame.Plot('Coupling') # can be Couplimg  -  STD
