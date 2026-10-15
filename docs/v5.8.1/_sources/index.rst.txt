PyMieSim documentation
======================

Calculate light scattering by spheres, core-shell spheres, and infinite
cylinders. Build a single-particle model, sweep experimental parameters, or
fit physical parameters to measured data with explicit units.

**Version:** |version|

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Get started
      :link: getting_started
      :link-type: doc

      Install PyMieSim and calculate scattering efficiency and collected power.

   .. grid-item-card:: Follow a workflow
      :link: workflows
      :link-type: doc

      Run parameter sweeps, fit measurements, average particle sizes, or plot fields.

   .. grid-item-card:: Explore examples
      :link: examples
      :link-type: doc

      Run complete scripts and learn to interpret their plots and numerical outputs.

   .. grid-item-card:: Look up the API
      :link: code/index
      :link-type: doc

      Find constructors, supported measures, and result types.

Choose your starting point
--------------------------

* New to the package? Start with :doc:`getting_started`.
* Have a calculation in mind? Use :ref:`choosing_an_interface`.
* Need to understand the physics? Read :doc:`theory` and :doc:`measures`.
* Installation failing? Check :doc:`compatibility` and :doc:`troubleshooting`.
* Preparing a reproducible study? See :doc:`reproducibility` and :doc:`performance`.

PyMieSim models isolated-particle scattering. Particle-size averaging assumes
independent scattering; it does not add interactions between particles.
See :doc:`workflows/particle_size_distributions` for the averaging assumptions.

.. toctree::
   :maxdepth: 2
   :hidden:

   Start <getting_started>
   Guide <guide>
   Examples <examples>
   API <code/index>
   Resources <resources>
