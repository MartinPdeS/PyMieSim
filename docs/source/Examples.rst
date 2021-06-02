
Examples
========


Creating Scatterers
-------------------

.. literalinclude:: ../../examples/Index.py
    :language: python
    :caption: **Scatterer: Index**
    :lines: 1-
    :linenos:



-------------------

.. literalinclude:: ../../examples/Material.py
    :language: python
    :caption: **Scatterer: Material**
    :lines: 1-
    :linenos:



-------------------------

.. literalinclude:: ../../examples/ScatProperties.py
    :language: python
    :caption: **Scatterer: Scatterer Properties**
    :lines: 1-
    :linenos:


-------------------------

.. literalinclude:: ../../examples/SamplesProperties.py
    :language: python
    :caption: **Scatterer: Samples properties**
    :lines: 1-
    :linenos:

----------------


Computing fields
----------------

.. literalinclude:: ../../examples/S1S2.py
    :language: python
    :linenos:
    :caption: **Scatterer: S1-S2**
    :lines: 1-


.. image:: ../images/S1S2.png
   :width: 600



-----------------


.. literalinclude:: ../../examples/Stokes.py
    :language: python
    :caption: **Scatterer: Stokes**
    :lines: 1-
    :linenos:

.. image:: ../images/Stokes.png
   :width: 600



-------------------------

.. literalinclude:: ../../examples/FarField.py
    :language: python
    :caption: **Scatterer: full far-field**
    :lines: 1-
    :linenos:

.. image:: ../images/FarField.png
   :width: 600



-------------------------

.. literalinclude:: ../../examples/SPF.py
    :language: python
    :caption: **Scatterer: phase function**
    :lines: 1-
    :linenos:



.. image:: ../images/SPF.png
   :width: 600


--------------------


Creating Detectors
------------------


.. literalinclude:: ../../examples/Photodiode.py
    :language: python
    :caption: **Detector: Photodiode**
    :lines: 1-
    :linenos:

.. image:: ../images/Photodiode.png
   :width: 600



----------------

.. literalinclude:: ../../examples/LPMode.py
    :language: python
    :caption: **Detector: LPMode**
    :lines: 1-
    :linenos:


.. image:: ../images/LPMode.png
   :width: 600


------------------------------


Computing coupling
------------------


.. literalinclude:: ../../examples/Scatterer-Photodiode.py
    :language: python
    :caption: **Coupling: Scatterer-Photodiode**
    :lines: 1-
    :linenos:


Output: (6.57e+01 nWatt)


--------------------------

Experiment: properties
----------------------


.. literalinclude:: ../../examples/Qscattering.py
    :language: python
    :caption: **ScattererSet: Qscattering**
    :lines: 1-
    :linenos:

.. image:: ../images/Qscattering.png
   :width: 600


------------------------------

.. literalinclude:: ../../examples/Qsca-vs-diameter.py
    :language: python
    :caption: **Experiment: Qsca-vs-diameter**
    :lines: 1-
    :linenos:

.. image:: ../images/Qsca-vs-diameter.png
  :width: 600


--------------------------

.. literalinclude:: ../../examples/Mie-resonances.py
    :language: python
    :caption: **Experiment: Mie-resonances**
    :lines: 1-
    :linenos:


.. image:: ../images/Mie-resonances.png
  :width: 600


---------------------------------


Experiment: coupling
--------------------

.. literalinclude:: ../../examples/Coupling-vs-diameter.py
    :language: python
    :caption: **Experiment: Coupling-vs-diameter**
    :lines: 1-
    :linenos:

.. image:: ../images/Coupling-vs-diameter.png
  :width: 600


-----------------------------------

.. literalinclude:: ../../examples/Coupling-vs-wavelength.py
    :language: python
    :caption: **Experiment: Coupling-vs-wavelength**
    :lines: 1-
    :linenos:

.. image:: ../images/Coupling-vs-wavelength.png
   :width: 600


------------------------------


Experiment: optimization
------------------------

.. literalinclude:: ../../examples/Opt-1-parameter.py
    :language: python
    :caption: **Optimization: Opt-1-parameter**
    :lines: 1-
    :linenos:


..
  .. image:: ../images/Opt1Param.png
     :width: 600


------------------------------

.. literalinclude:: ../../examples/Opt-2-parameter.py
    :language: python
    :caption: **Optimization: Opt-2-parameter**
    :lines: 1-
    :linenos:

------------------------------

=========================================================
Plot the scattering efficiency of a sphere using PyMieSim
=========================================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter. The graph is
the following:

.. image:: ../images/ScatteringEfficiency.png

Making this graph using PyMieSim is very simple. Since the graph is made with
matplotlib and numpy, the first step is to make sure these packages are installed
and then to import them:

.. literalinclude:: ../../examples/ScatteringEfficiency.py
    :language: python
    :caption: **Import librairies**
    :lines: 2-11
    :linenos:

The following step is to create a light source and an empty list that will
contain the heatmap values:

.. literalinclude:: ../../examples/ScatteringEfficiency.py
    :language: python
    :caption: **Add light source**
    :lines: 13-18
    :linenos:

The next step is the most important since in generates the list of lists that
will be used to make the heatmap. It consists of two for loops that loop through
different values of Diameter and Index which are two parameters of a Sphere in
PyMieSim. Note that some scale adjustments were made to adjust the graph's aspect
ratio:

.. literalinclude:: ../../examples/ScatteringEfficiency.py
    :language: python
    :caption: **Create list of lists**
    :lines: 20-38
    :linenos:

Following this step, we have a list of lists named heatmap that will generate the
graph. To do so, it is first converted in a numpy array and then some operations
are made using matplotlib to plot the graph accurately:

.. literalinclude:: ../../examples/ScatteringEfficiency.py
    :language: python
    :caption: **Plot the graph**
    :lines: 40-69
    :linenos:

The entire code to generate the graph of the scattering efficiency of a sphere
can be found below in a single block. The result should be the plot shown under
the code block.

.. literalinclude:: ../../examples/ScatteringEfficiency.py
    :language: python
    :caption: **Scattering efficiency graph**
    :lines: 1-
    :linenos:

.. image:: ../images/ScatteringEfficiency.png
