.. index::
  single: Green's functions; block Green's function
  module: GFBloc_ImFreq

.. _GFBloc_ImFreq:

Matsubara Green's function in imaginary frequencies (GFBloc_ImFreq)
=====================================================================================

This object stores a matrix valued Green function in imaginary frequencies
in the Matsubara formalism : 

.. math::
  G_{\alpha \beta} (i \omega_n) \equiv \int_0^\beta G_{\alpha \beta} ( \tau) e^{-i \omega_n \tau}

where :math:`G(\tau)` is the :ref:`Matsubara Green function in imaginary time <GFBloc_ImTime>`.


Reference
---------------

.. autoclass:: pytriqs.base.GF_Local.GFBloc_ImFreq
  :members: density, setFromFourierOf, setFromLegendre
 

HDF5 data scheme
^^^^^^^^^^^^^^^^^^^^^^

The GFBloc_ImFreq (TRIQS_HDF5_data_scheme= "GFBloc_ImFreq") is decomposed in the following objects : 

=========================   ===========================  ===========================================================================
Name                        Type                         Meaning
=========================   ===========================  ===========================================================================
Mesh                        MeshGF                       The mesh
Tail                        TailGF                       The tail
Data                        3d numpy of complex          Data[i1,i2,n] is the element of the Green function where :
                                                          * i1, i2 are the indices
                                                          * n is the Matsubara index 
                                                         
                                                         For real GF in time, we store only starting from n=0 and the
                                                         rest is deduced by symmetry
IndicesL,IndicesR           string                       The Python repr of the indices, e.g. (1,2), or (1,)
                                                         repr(this_string) reproduces the indices 
Name                        string                       Name of the Green function block
Note                        string                       Note 
=========================   ===========================  ===========================================================================


Example 
--------------- 

.. plot:: green/block/green_imfreq.py
   :include-source:
   :scale: 70 
