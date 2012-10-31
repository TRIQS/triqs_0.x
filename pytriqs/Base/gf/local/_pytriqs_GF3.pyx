#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
#from libcpp.vector cimport vector
cimport cython  
from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
  
# -------------------- Some generic tools -------------------------------
 
cdef extern from "triqs/gf/tools.hpp" namespace "triqs::gf" : 
    cdef enum statistic_enum "triqs::gf::statistic_enum" :
        Boson,Fermion 
    
    cdef cppclass nothing : 
        nothing ()

    cdef cppclass freq_infty: 
        freq_infty ()

    cdef cppclass mesh_pt_generator "triqs::gf::mesh_pt_generator" [MeshType] :  
        mesh_pt_generator( MeshType * )
        mesh_pt_generator()
        complex to_point()
        mesh_pt_generator operator++()
        bint at_end()
        void increment() 
 
# ----------- Local high frequency tail --------------------------

cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail_view "triqs::gf::local::tail_view"  :
        tail_view()
        tail_view(array_view[dcomplex,THREE,COrder] , int) #except +
        void rebind (tail_view &)

cdef class TailGF_c:
    cdef tail_view _c
    def __init__(self, a, int omin):
        self._c.rebind(  tail_view( array_view[dcomplex,THREE,COrder](a), omin))

include "./matsubara_freq.pyx"
include "./matsubara_time.pyx"  
