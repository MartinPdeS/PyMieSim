Scatterer material
==================

Built-in material
-----------------

PyMieSim allows you to use the bank of **Material** to emulate scatterer refracive index which vary
as a function of the source wavelength.

To use the **Material** follow the snippet:

.. code-block:: python
  :linenos:

  from PyMieSim.Scatterer import Sphere
  from PyMieSim.Source import PlaneWave
  from PyMieSim import Material

  Source = PlaneWave(Wavelength   = 450e-9,
                    Polarization = 0,
                    E0           = 1)

  Scat = Sphere(Diameter    = 800e-9,
               Source       = Source,
               Material     = Material('BK7'))

----

Generate new material
----------------------

PyMieSim has a toolbox to fetch and load new material data.
Using the largest refracitve index bank: `refractiveindex.info <https://refractiveindex.info/>`_  PyMieSim can directly
download and save complex refractive index for a specific material.
Doing so the toolbox will save locally the data in the **PyMieSim/Data/_Material/data** folder and will update the
**PyMieSim/Meta.json** file.
Here is what the Meta file look like


.. literalinclude:: ../../PyMieSim/Data/_Material/Meta.json
    :language: JSON
    :linenos:


----

 In order to generate new Material to use with PyMieSim one can use the following snippet



 .. literalinclude:: ../../examples/Extra:New-Material-BK7.py
     :language: python
     :linenos:


.. image:: ../images/Extra:New-Material-BK7.png
   :width: 600


After executing this code a new entry to **Meta.json** will be added or updated.
Afterward the user can use this new material for my computation.


Here is another examples



.. literalinclude:: ../../examples/Extra:New-Material-Silver.py
    :language: python
    :linenos:


.. image:: ../images/Extra:New-Material-Silver.png
  :width: 600
