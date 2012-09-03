The Monte Carlo class
=====================

The TRIQS library has a class called ``mc_generic`` which allows you to write
Monte Carlo algorithm in a simple framework. The class takes care of the basic
mechanics which is common to any Monte Carlo method so that you can focus on
the implementation details of your specific algorithm.


A generic Monte Carlo
---------------------

Here is a diagram of a generic Monte Carlo:

.. image:: loop.png
   :width: 700
   :align: center


..
  * Easily add new moves and measures
  * For complex MC algorithms, with 6 or 7 different moves, it helps to clearly separate the code for each move.
  * Parallelism is automatic.
  * The random generator is a simple parameter, it can be chosen dynamically. 

..
  The `mc_generic` class is a generic version of the algorithms, with `moves` and `measures`.
  The user  : 
    
    - writes move classes, modelling the Move concept.
    - writes measure classes, modelling the Measure concepts.
    - register them in a mc_generic object MC
    - call `MC.run(...)`  ... and that is (almost) it  !

