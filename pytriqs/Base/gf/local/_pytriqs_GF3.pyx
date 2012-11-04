#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
cimport cython  
from gf_matsubara_freq cimport *
from gf_matsubara_time cimport *
from fourier_matsubara cimport *
from libcpp.vector cimport vector

cdef class TailGF_c:
    cdef tail_view_c _c
    def __init__(self, a, int omin):
        self._c =  tail_view_c( array_view[dcomplex,THREE,COrder](a), omin)

include "./matsubara_freq.pyx"
include "./matsubara_time.pyx"  


cdef extern from "pytriqs/Base/gf/local/a.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_view_c "triqs::gf::gf_view<triqs::gf::block<matsubara_freq> >" :
        gf_block_view_c()
    
    cdef gf_block_view_c make_gf_block_view_c "triqs::gf::block<matsubara_freq>::make_gf_view" (  vector[gf_view_freq_c] &) 

    cdef void test_block_c "test_block" ( gf_block_view_c&)  

def test_block (G) : 
    cdef vector[gf_view_freq_c] v
    for item in G:
        v.push_back((<GFBloc_ImFreq_cython>item)._c)
    
    cdef gf_block_view_c GG = make_gf_block_view_c (v)  
    test_block_c (GG)
    
    test_block_c (make_gf_block_view_c (v) )




