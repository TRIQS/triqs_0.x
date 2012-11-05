#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
cimport cython  
from gf_matsubara_freq cimport *
from gf_matsubara_time cimport *
from fourier_matsubara cimport *
from arrays cimport *   
from libcpp.vector cimport vector
#from tail cimport *

cdef class TailGF_c:
    #cdef tail_view_c _c
    def __init__(self, a, int omin):
        self._c =  tail_view_c( array_view[dcomplex,THREE,COrder](a), omin)

include "./matsubara_freq.pyx"
include "./matsubara_time.pyx"  

import numpy

cdef extern from "pytriqs/Base/gf/local/a.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_view_c "triqs::gf::gf_view<triqs::gf::block<matsubara_freq> >" :
        gf_block_view_c()
    
    cdef gf_block_view_c make_gf_block_view_c "triqs::gf::block<matsubara_freq>::make_gf_view" (  vector[gf_view_freq_c] &) 

    cdef void test_gf_c "test_gf" ( gf_view_freq_c &)  
    cdef void test_block_c "test_block" ( gf_block_view_c&)  
    cdef void test2_c "test2" ( matrix_view[double,COrder] &) except + 
    cdef void test3_c "test3" ( vector[matrix_view[double,COrder] ] &) except + 
    #cdef void test3_c "test3" ( vector[double] &)  

def test_block (G) : 
    #cdef vector[gf_view_freq_c] v_c
    #for item in G:
    #    #v_c.push_back((<GFBloc_ImFreq_cython>item)._c)
    #    v_c.push_back(to_C_GFBloc_ImFreq(item))
    
    #cdef gf_block_view_c GG = make_gf_block_view_c (v_c)  
    test_block_c (to_gf_block_view_c(G)) 
    
     
    #test_gf_c (G[0]) 
def test2():
    a = numpy.array([[1,2],[3,4]])
    test2_c(a)

def test3():
    a = numpy.array([[1.0,2],[3,4]])
    w = [a,a]
    test3_c(w)
    print a


