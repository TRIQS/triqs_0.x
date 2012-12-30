
.. _hilbert_transform:

.. module:: pytriqs.base.dos.hilbert_transform


Hilbert Transform
=======================================
TRIQS comes with a Hilbert transform. Let us look at an example:

.. literalinclude:: ex_Hilbert.py


Given a density of states `d` (here for a tight-binding model), the Hilbert transform `H` is defined is defined in the following way::

    H = Hilbert_Transform(d)

To construct a Green's function::

    G = GFBloc_ImFreq(Indices = ['up','down'], Beta = 20)
    Sigma0 = GFBloc_ImFreq(Indices = ['up','down'], Beta = 20); Sigma0.zero()
    G <<= H(Sigma = Sigma0,mu=0.)


.. autoclass:: pytriqs.base.dos.Hilbert_Transform
   :members: __call__
   :undoc-members:
  

