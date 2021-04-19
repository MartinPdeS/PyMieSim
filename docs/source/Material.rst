
PyMieSim has a toolbox to fetch and load new material data.
Using the largest refracitve index bank: `refractiveindex.info <https://refractiveindex.info/>`_  PyMieSim can directly
download and save complex refractive index for a specific material.
Doing so the toolbox will save locally the data in the **PyMieSim/_Material/data** folder and will update the
**PyMieSim/Meta.json** file.
Here is what the Meta file look like


.. literalinclude:: ../../PyMieSim/_Material/Meta.json
    :language: JSON
    :linenos:


----

 In order to generate new Material to use with PyMieSim one can use the following snippet



 .. literalinclude:: ../../tests/Examples/New-Material.py
     :language: python
     :linenos:


.. image:: ../images/New-Material.png
   :width: 600


After executing this code a new entry to **Meta.json** will be added or updated.
Afterward the user can use this new material for my computation.


Here is another examples



.. literalinclude:: ../../tests/Examples/New-Material-Silver.py
    :language: python
    :linenos:


.. image:: ../images/New-Material-Silver.png
  :width: 600
