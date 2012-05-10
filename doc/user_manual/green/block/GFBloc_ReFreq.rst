.. index::
  single: Green's functions; block Green's function
  module: GFBloc_ReFreq

.. _GFBloc_ReFreq:

Retarded Green's function in real frequencies (GFBloc_ReFreq)
=====================================================================================

This object stores a matrix valued Green function in real frequencies : 

.. math::
  G_{\alpha \beta} (\omega) \equiv \int_{-\infty}^{\infty} G_{\alpha \beta} ( t ) e^{-i \omega t}

where :math:`G(t)` is the :ref:`Green function in real time <GFBloc_ReTime>`.


Reference
---------------

.. autoclass:: pytriqs.Base.GF_Local.GFBloc_ReFreq
  :members: density, setFromFourierOf, transpose, conjugate, setFromPadeOf
 


HDF5 data scheme
^^^^^^^^^^^^^^^^^^^^^^

The GFBloc_ReFreq (TRIQS_HDF5_data_scheme = "GFBloc_ReFreq") is decomposed in the following objects : 

=========================   ===========================  ===========================================================================
Name                        Type                         Meaning
=========================   ===========================  ===========================================================================
Mesh                        MeshGF                       The mesh
Tail                        TailGF                       The tail
Data                        3d numpy of complex          Data[i1,i2,n] is the element of the Green function where :
                                                          * i1, i2 are the indices
                                                          * n is the frequency index 
IndicesL,IndicesR           string                       The Python repr of the indices, e.g. (1,2), or (1,)
                                                         repr(this_string) reproduces the indices 
Name                        string                       Name of the Green function block
Note                        string                       Note 
=========================   ===========================  ===========================================================================


Examples 
---------------

.. plot:: green/block/green_refreq.py
   :include-source:
   :scale: 70

Note that `g` is a **retarded** Green's function.


The next example demonstrates how a real frequency Green's function can be
reconstructed from an imaginary frequency counterpart using setFromPadeOf()
method. Note that in order to use this method you need to activate the
``Use_Pade`` option when you run cmake (``cmake -DUse_Pade=ON``).

.. plot:: green/block/green_pade.py
   :include-source:
   :scale: 70
