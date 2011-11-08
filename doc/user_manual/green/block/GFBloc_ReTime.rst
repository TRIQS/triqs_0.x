.. index::
  single: Green's functions; block Green's function
  module: GFBloc_ReTime

.. _GFBloc_ReTime:

Retarded Green's function in real time (GFBloc_ReTime)
=====================================================================================

This object stores a retarded matrix valued Green function in real time : 

.. math::
  G_{\alpha \beta} (t) \equiv - \theta(t) <\{c_{\alpha}(t) c^{\dagger}_{\beta} (0)\}> 


Reference
---------------

.. autoclass:: pytriqs.Base.GF_Local.GFBloc_ReTime
  :members: density, setFromInverseFourierOf, transpose, conjugate
 

HDF5 data scheme
^^^^^^^^^^^^^^^^^^^^^^

The GFBloc_ReTime (TRIQS_HDF5_data_scheme = "GFBloc_ReTime") is decomposed in the following objects : 

=========================   ===========================  ===========================================================================
Name                        Type                         Meaning
=========================   ===========================  ===========================================================================
Mesh                        MeshGF                       The mesh
Tail                        TailGF                       The tail
Data                        3d numpy of complex          Data[i1,i2,n] is the element of the Green function where :
                                                          * i1, i2 are the indices
                                                          * n is the index of the time slice
IndicesL,IndicesR           string                       The Python repr of the indices, e.g. (1,2), or (1,)
                                                         repr(this_string) reproduces the indices 
Name                        string                       Name of the Green function block
Note                        string                       Note 
=========================   ===========================  ===========================================================================


Example 
---------------

.. plot:: green/block/green_retime.py
   :include-source:
   :scale: 70
  
