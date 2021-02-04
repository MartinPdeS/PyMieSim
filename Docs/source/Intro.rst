Introduction
============

This project (SuPyModes) aims to create a python library that gives the necesary tool to simulate propagation mode of light inside fibers coupler.


Getting started
---------------



Here is the list of packages you need to use this library:
    - Python (3+)
    - Numpy
    - Scipy
    - Pandas
    - Fibermodes: https://github.com/cbrunet/fibermodes
    - Mayavi


Mie Theory
----------

Calculates S\ :sub:`1` and S\ :sub:`2` at μ=cos(θ), where θ is the scattering angle.

 S\ :sub:`1` and S\ :sub:`2` are calculated by:

             :math:`{\displaystyle S_1=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\pi_n+b_n\tau_n)}`

             :math:`{\displaystyle S_2=\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}(a_n\tau_n+b_n\pi_n)}`


 Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses :py:func:`Mie_ab` to calculate :math:`a_n` and :math:`b_n`, and then calculates *Q* via:

              :math:`${\displaystyle Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)}$`





Run test simulation
-------------------

Here an example of command to run a simulation on linux (ubuntu 16.04).

.. code-block:: console
   :linenos:
