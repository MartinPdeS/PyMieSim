Examples
========

Scatterer: S1-S2
----------------

.. code-block:: console
   :linenos:

   from PyMieSim.Scatterer import Sphere
   from PyMieSim.Source import PlaneWave

   LightSource = PlaneWave(Wavelength = 450e-9,
                           Polarization = 0)

   Scat = Sphere(Diameter    = 300e-9,
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

   from PyMieSim.Scatterer import Sphere
   from PyMieSim.Source import PlaneWave

   LightSource = PlaneWave(Wavelength = 450e-9,
                           Polarization = 0)

   Scat = Sphere(Diameter    = 300e-9,
                 Source      = LightSource,
                 Index       = 1.4)


   Fields = Scat.FarField(Num=100)

   Fields.Plot()


.. image:: ../images/Fields.png
   :width: 600

Scatterer: phase function
-------------------------

.. code-block:: console
   :linenos:

   from PyMieSim.Scatterer import Sphere
   from PyMieSim.Source import PlaneWave

   LightSource = PlaneWave(Wavelength = 450e-9,
                          Polarization = 0)

   Scat = Sphere(Diameter    = 800e-9,
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

   from PyMieSim.Source import PlaneWave
   from PyMieSim.Detector import Photodiode

   LightSource = PlaneWave(Wavelength = 450e-9,
                           Polarization = 0)

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

   from PyMieSim.Source import PlaneWave
   from PyMieSim.classes.Detector import LPmode

   LightSource = PlaneWave(Wavelength = 450e-9,
                           Polarization = 0)

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

   from PyMieSim.Source import PlaneWave
   from PyMieSim.Detector import LPmode
   from PyMieSim.Scatterer import Sphere

   LightSource = PlaneWave(Wavelength = 450e-9,
                           Polarization = 0)

   Detector = LPmode(Mode         = (0, 1,'h'),
                     Sampling     = 201,
                     NA           = 0.2,
                     GammaOffset  = 0,
                     PhiOffset    = 0,
                     CouplingMode = 'Centered')


   Scat = Sphere(Diameter    = 300e-9,
                 Source      = LightSource,
                 Index       = 1.4)

   Coupling = Detector.Coupling(Scatterer = Scat)

   print(Coupling)


Output: (6.70121870391961e-18)


ScattererSet: Qscattering
--------------------------

.. code-block:: console
   :linenos:

   import numpy as np
   from PyMieSim.Source import PlaneWave
   from PyMieSim.Sets import ScattererSet

   LightSource = PlaneWave(Wavelength = 450e-9,
                          Polarization = 0)


   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 15000e-9, 400),
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
   from PyMieSim.Source import PlaneWave
   from PyMieSim.Detector import LPmode
   from PyMieSim.Sets import ScattererSet, ExperimentalSet

   LightSource = PlaneWave(Wavelength = 450e-9, Polarization = 0)



   Detector0 = LPmode(NA               = 0.2,
                     Sampling          = 401,
                     GammaOffset       = 0,
                     PhiOffset         = 20,
                     Mode              = (0,1),
                     CouplingMode      = 'Mean')

   Detector1 = LPmode(NA               = 0.2,
                     Sampling          = 401,
                     GammaOffset       = 0,
                     PhiOffset         = 20,
                     Mode              = (1,1),
                     CouplingMode      = 'Mean')





   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 1500e-9, 500),
                          RIList        = np.linspace(1.5, 1.5, 1).round(1),
                          Source        = LightSource)





   Set = ExperimentalSet(ScattererSet  = ScatSet,
                         Detectors     = (Detector0, Detector1))


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
  from PyMieSim.Detector import Photodiode, LPmode
  from PyMieSim.Source import PlaneWave
  from PyMieSim.Optimizer import Simulator
  from PyMieSim.Sets import ExperimentalSet, ScattererSet

  LightSource = PlaneWave(Wavelength = 450e-9,
                         Polarization = 0)

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
