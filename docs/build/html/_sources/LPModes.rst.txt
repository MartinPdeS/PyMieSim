
PyMieSim offer to the user to create phase-sensitive detectors such as fiber optics LP modes.
To check how to use them, please refer to the examples section

A few LP mode are already precomputed but it might be convenient to add mores.
To do so read the following.

LP-modes: Generate files
------------------------

I have prealably compilated 12 LP-modes which are:

1. LP01
2. LP11
3. LP21
4. LP02
5. LP31
6. LP12
7. LP41
8. LP22
9. LP03
10. LP51
11. LP32
12. LP13

So you can already use them, no need to reproduce. However if you want to
use another mode you first need to install the fibermodes package using the
following command:

.. code-block:: python
  :linenos:

  pip install https://github.com/cbrunet/fibermodes.git

Then you can use the PyMieSIm FiberModes module as follow:

.. code-block:: python
  :linenos:

  from PyMieSim.FiberModes import Genfiles

  Genfiles([(5,2)], padWidth = 2000, Num = 251)


Here **(5,2)** is the mode we want to add, **padWidth** is the pre-FFT zero-padding,
it controls the width of the computed LP mode it should stay around 2000.
Finally **Num** is the lateral resolution for the LP mode. It is safe to keep it at 251.
Right now **Num** has to be an odd number. Maybe I will look for it soon but not now, sorry! 
