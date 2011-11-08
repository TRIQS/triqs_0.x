.. _random:


Random number generators
=========================

It is possible (and highly recommended!) to use different random number generators with the CTQMC solver.
A list of the available generators is obtained with the following script::

  from pytriqs.Solvers.HybridizationExpansion import Solver
  print Solver.Random_Generators_Available()

The names in this list can then be used to set the ``Random_Generator_Name`` keyword argument
in the constructor of the CTQMC solver. Here is an example where the same run is done
twice, but using different random generator.

.. literalinclude:: several_random.py

The resulting Green's functions are:
  
.. image:: random.png
   :width: 700
   :align: center

As one can see from the plots, the two random number generators lead to very similar results.
On these short runs, one also notices that more statistics will be needed to have a reliable
Green's function at low frequencies.
