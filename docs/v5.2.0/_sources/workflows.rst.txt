.. _workflows:

Workflows
=========

This section is organized around common tasks rather than individual classes.
Most examples use the public top-level API, so the same patterns can be used
without importing generated C++ extension modules.

.. toctree::
   :maxdepth: 1

   workflows/first_simulation
   workflows/parameter_sweeps
   workflows/detector_coupling
   workflows/visualizing_fields

The typical PyMieSim workflow is:

.. code-block:: text

   Source + material → scatterer → Simulation / Experiment → result

For runnable examples, see the :ref:`examples_gallery`.
