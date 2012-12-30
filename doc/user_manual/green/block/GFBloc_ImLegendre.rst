.. index::
  single: Green's functions; block Green's function
  module: GFBloc_ImLegendre

.. _GFBloc_ImLegendre:

Matsubara Green's function with Legendre polynomials (GFBloc_ImLegendre)
=====================================================================================

This object stores a matrix-valued Matsubara Green's function in imaginary time using
Legendre polynomials. It is connected to the usual imaginary-time Matsubara Green's
function through:

.. math::
  G_{l, \alpha \beta} = \sqrt{2l+1} \int_0^\beta d\tau P_l (x(\tau)) G_{\alpha \beta}(\tau),


where :math:`x(\tau) = 2\tau/\beta -1` and the :math:`P_l` are the Legendre
polynomials.  :math:`G_{\alpha \beta}(\tau)` is the :ref:`Matsubara Green's
function in imaginary time <GFBloc_ImTime>`.


Reference
---------------

.. autoclass:: pytriqs.base.GF_Local.GFBloc_ImLegendre
  :members: density, determine_tail, enforce_discontinuity, copy_and_truncate
 

HDF5 data scheme
^^^^^^^^^^^^^^^^^^

The GFBloc_ImLegendre (TRIQS_HDF5_data_scheme= "GFBloc_ImLegendre") is decomposed in the following objects : 

=========================   ===========================  ===========================================================================
Name                        Type                         Meaning
=========================   ===========================  ===========================================================================
Mesh                        MeshGF                       The mesh
Tail                        TailGF                       The tail
Data                        3d numpy of complex          Data[i1,i2,n] is the element of the Green function where :
                                                          * i1, i2 are the indices
                                                          * n is the Legendre coefficient index 
IndicesL, IndicesR          string                       The Python repr of the indices, e.g. (1,2), or (1,)
                                                         repr(this_string) reproduces the indices 
Name                        string                       Name of the Green's function block
Note                        string                       Note
=========================   ===========================  ===========================================================================

Example
--------------- 

.. plot:: green/block/green_imlegendre.py
   :include-source:
   :scale: 70 

