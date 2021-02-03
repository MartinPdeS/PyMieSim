Examples
========

Scatterer: S1-S2
----------------

.. code-block:: console
   :linenos:
   
   from PyMieCoupling.classes.Scattering import Scatterer
   from PyMieCoupling.utils import Source
   
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
   from PyMieCoupling.utils import Source
   
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
   from PyMieCoupling.utils import Source
   
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
   
   from PyMieCoupling.utils import Source
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
   
   from PyMieCoupling.utils import Source
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
   
   
   
Coupling: Scatterer-LPMode
--------------------------

.. code-block:: console
   :linenos:
   
   from PyMieCoupling.utils import Source
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
                    
                    
   Coupling = Detector.Coupling(Scatterer = Scat)  
   
   print(Coupling)
   
   
   
ScattererSet: Qscattering
--------------------------

.. code-block:: console
   :linenos:
   
   import numpy as np
   from PyMieCoupling.utils import Source
   from PyMieCoupling.classes.Sets import ScattererSet

   LightSource = Source(Wavelength   = 950e-9,
                        Polarization = 0)


   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 500e-9, 100),
                          RIList        = np.linspace(1.5, 1.5, 1).round(1),
                          Source        = LightSource) 
   
     
   Qsca = ScatSet.Qsca()
   
   Qsca.Plot()
   
   
   
   
ExperimentalSet: Coupling
----------------------------

.. code-block:: console
   :linenos:
   
   import numpy as np
   from PyMieCoupling.utils import Source
   from PyMieCoupling.classes.Sets import ScattererSet

   LightSource = Source(Wavelength   = 950e-9,
                        Polarization = 0)


   ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 500e-9, 100),
                          RIList        = np.linspace(1.5, 1.5, 1).round(1),
                          Source        = LightSource) 
   
     
   Qsca = ScatSet.Qsca()
   
   Qsca.Plot()
      
   
   
   
   
   
   
   
   
   
   
   
