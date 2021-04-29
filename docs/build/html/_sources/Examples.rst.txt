
Examples
========


Scatterer: Index
-------------------

.. literalinclude:: ../../examples/Index.py
    :language: python
    :linenos:


Scatterer: Material
-------------------

.. literalinclude:: ../../examples/Material.py
    :language: python
    :linenos:



Scatterer: S1-S2
----------------




.. literalinclude:: ../../examples/S1S2.py
    :language: python
    :linenos:





.. image:: ../images/S1S2.png
   :width: 600


-----



Scatterer: Stokes
-----------------

.. literalinclude:: ../../examples/Stokes.py
    :language: python
    :linenos:

.. image:: ../images/Stokes.png
   :width: 600

-----



Scatterer: full far-field
-------------------------

.. literalinclude:: ../../examples/FarField.py
    :language: python
    :linenos:

.. image:: ../images/FarField.png
   :width: 600

-----



Scatterer: phase function
-------------------------

.. literalinclude:: ../../examples/SPF.py
    :language: python
    :linenos:



.. image:: ../images/SPF.png
   :width: 600

-----

Detector: Photodiode
--------------------

.. literalinclude:: ../../examples/Photodiode.py
    :language: python
    :linenos:

.. image:: ../images/Photodiode.png
   :width: 600

-----




Detector: LPMode
----------------

.. literalinclude:: ../../examples/LPMode.py
    :language: python
    :linenos:


.. image:: ../images/LPMode.png
   :width: 600

-----




Coupling: Scatterer-Photodiode
------------------------------

.. literalinclude:: ../../examples/Scatterer-Photodiode.py
    :language: python
    :linenos:


Output: (6.57e+01 nWatt)

-----




ScattererSet: Qscattering
--------------------------

.. literalinclude:: ../../examples/Qscattering.py
    :language: python
    :linenos:

.. image:: ../images/Qscattering.png
   :width: 600


-----




Experiment: Qsca-vs-diameter
------------------------------

.. literalinclude:: ../../examples/Qsca-vs-diameter.py
    :language: python
    :linenos:

.. image:: ../images/Qsca-vs-diameter.png
  :width: 600

-----




Experiment: Mie-resonances
--------------------------

.. literalinclude:: ../../examples/Mie-resonances.py
    :language: python
    :linenos:


.. image:: ../images/Mie-resonances.png
  :width: 600

-----




Experiment: Coupling-vs-diameter
---------------------------------

.. literalinclude:: ../../examples/Coupling-vs-diameter.py
    :language: python
    :linenos:

.. image:: ../images/Coupling-vs-diameter.png
  :width: 600


-----

Experiment: Coupling-vs-wavelength
-----------------------------------

.. literalinclude:: ../../examples/Coupling-vs-wavelength.py
    :language: python
    :linenos:

.. image:: ../images/Coupling-vs-wavelength.png
   :width: 600

-----

Optimization: Opt-1-parameter
------------------------------

.. literalinclude:: ../../examples/Opt-1-parameter.py
    :language: python
    :linenos:


**Output:**


| Call Number : 1             	 PhiOffset: 1.00000e-01             	 Result: -7.3947105131e-03
| Call Number : 2             	 PhiOffset: 1.01000e+01             	 Result: -4.5216666010e-03
| Call Number : 3             	 PhiOffset: -9.90000e+00             	 Result: -4.6103038869e-03
| Call Number : 4             	 PhiOffset: -4.90000e+00             	 Result: -6.5239220916e-03
| Call Number : 5             	 PhiOffset: 2.60000e+00             	 Result: -7.1347913938e-03
| Call Number : 6             	 PhiOffset: -1.15000e+00             	 Result: -7.3444635289e-03
| Call Number : 7             	 PhiOffset: 7.25000e-01             	 Result: -7.3742571154e-03
| Call Number : 8             	 PhiOffset: -2.12500e-01             	 Result: -7.3935287699e-03
| Call Number : 9             	 PhiOffset: 2.56250e-01             	 Result: -7.3924460072e-03
| Call Number : 10             	 PhiOffset: 2.18750e-02             	 Result: -7.3951290146e-03
| Call Number : 11             	 PhiOffset: -5.62500e-02             	 Result: -7.3950715369e-03
| Call Number : 12             	 PhiOffset: -1.71875e-02             	 Result: -7.3951597748e-03
| Call Number : 13             	 PhiOffset: -5.62500e-02             	 Result: -7.3950715369e-03
| Call Number : 14             	 PhiOffset: -3.67188e-02             	 Result: -7.3951305305e-03
| Call Number : 15             	 PhiOffset: -7.42188e-03             	 Result: -7.3951632409e-03
| Call Number : 16             	 PhiOffset: 2.34375e-03             	 Result: -7.3951592695e-03
| Call Number : 17             	 PhiOffset: -2.53906e-03             	 Result: -7.3951621849e-03
| Call Number : 18             	 PhiOffset: -9.86328e-03             	 Result: -7.3951630716e-03
| Call Number : 19             	 PhiOffset: -6.20117e-03             	 Result: -7.3951631512e-03
| Call Number : 20             	 PhiOffset: -8.03223e-03             	 Result: -7.3951632421e-03
| Call Number : 21             	 PhiOffset: -8.33740e-03             	 Result: -7.3951632319e-03
| Call Number : 22             	 PhiOffset: -7.87964e-03             	 Result: -7.3951632446e-03
| Call Number : 23             	 PhiOffset: -7.72705e-03             	 Result: -7.3951632451e-03
| Call Number : 24             	 PhiOffset: -7.57446e-03             	 Result: -7.3951632439e-03
| Call Number : 25             	 PhiOffset: -7.82705e-03             	 Result: -7.3951632450e-03
| fun: -0.007395163244966126
| maxcv: 0.0
| message: 'Optimization terminated successfully.'
| nfev: 25
| status: 1
| success: True
| x: array([-0.00782705])


..
  .. image:: ../images/Opt1Param.png
     :width: 600


-----

Optimization: Opt-2-parameter
------------------------------

.. literalinclude:: ../../examples/Opt-2-parameter.py
    :language: python
    :linenos:
