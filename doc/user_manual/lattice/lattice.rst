.. index:: lattices

.. module:: pytriqs.lattice

.. _lattice:

Lattice Tools
=============

TRIQS comes with a couple of modules that allow to work with free fermions on Bravais lattices, with 
tight-binding hoppings.

.. warning::
   This part of the documentation largely is work in progress ...

Introductory example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ex1.py
   :lines: 1-34

Here, we take the simple square lattice in d=2 and compute the density of state. It produces : 

.. image:: ex1.png
   :width: 700
   :align: center

Reference manual
^^^^^^^^^^^^^^^^^^^^^^
 

.. toctree::

  bravais
  tightbinding
  dos
  hilbert
  sumk


