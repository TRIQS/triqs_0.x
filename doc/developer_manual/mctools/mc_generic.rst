
The mc_generic class
============================

The MonteCarlo class runs a Metropolis algorithm, given :

  - a set of Move with their proposition probabilities.
  - a set of Measure to handle measurements.


Algorithms
-------------------------


The algorithm runs by the class is::

  while Number_of_Cycles < Number_of_Cycles_max and  time < timelimit: 
    
   repeat Length_MC_Cycle times :  [ the Cycle]
     - choose a Move with its proposition probability
     - probability to accept = Move->Try()
     - Metropolis choice : signe = move_sign * accept_sign where
     - move_sign is the sign of Move->Try() and
     - accept_sign is the sign of Move->Accept()
   for all Measures : accumulate
  for all Measures : Finalize


Construction
-------------------------



 * The parameter dictionnary can be of alps::mcparams type or a python dictionnary

Doxygen documentation
--------------------------

The :mctools_doxy:`full C++ documentation<triqs::mc_tools::mc_generic>` is available here.

Breathe Documentation 
----------------------------------

.. doxygenclass:: triqs::mc_tools::mc_generic
  :project: mc_tools
  :members:
   
