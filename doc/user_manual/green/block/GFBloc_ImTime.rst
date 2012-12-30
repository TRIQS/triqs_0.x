.. index::
  single: Green's functions; block Green's function
  module: GFBloc_ImTime

.. _GFBloc_ImTime:

Matsubara Green's function in imaginary time (GFBloc_ImTime)
==============================================================================

This object stores a matrix valued Green function in imaginary time
in the Matsubara formalism : 

.. math::

  G_{\alpha \beta} (\tau) \equiv - \mathcal{T} <c_{\alpha}(\tau) c^{\dagger}_{\beta} (0)> \qquad \text{ for } 0 \leq \tau \leq \beta

 
Reference
---------------

.. autoclass:: pytriqs.base.GF_Local.GFBloc_ImTime
  :members: setFromInverseFourierOf, setFromLegendre, transpose, conjugate


HDF5 data scheme
^^^^^^^^^^^^^^^^^^^^^^

The GFBloc_ImTime (TRIQS_HDF5_data_scheme = "GFBloc_ImTime") is decomposed in the following objects : 

=========================   ===========================  ===========================================================================
Name                        Type                         Meaning
=========================   ===========================  ===========================================================================
Mesh                        MeshGF                       The mesh
Tail                        TailGF                       The tail
Data                        3d numpy of double           Data[i1,i2,n] is the element of the Green function where :
                                                          * i1, i2 are the indices
                                                          * n is the index of the time slice, starting at tau =0
IndicesL,IndicesR           string                       The Python repr of the indices, e.g. (1,2), or (1,)
                                                         repr(this_string) reproduces the indices 
Name                        string                       Name of the Green function block
Note                        string                       Note 
=========================   ===========================  ===========================================================================

Example
-------

.. plot:: green/block/green_imtime.py
   :include-source:
   :scale: 70 

