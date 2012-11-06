#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
cimport cython  

###############  Fourier  #########################

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_im_freq_c lazy_fourier          (gf_im_time_c & )
    gf_im_time_c lazy_inverse_fourier  (gf_im_freq_c & )

###############  Implementation de Tail  #########################

cdef class TailGF_c:
    #cdef tail_view_c _c
    def __init__(self, a, int omin):
        self._c =  tail_view_c( array_view[dcomplex,THREE,COrder](a), omin)

include "./matsubara_freq.pyx"
include "./matsubara_time.pyx"  

#---------------------- TEST FOR DEBUG ONLY ...------------------------------------------

cdef extern from "pytriqs/Base/gf/local/a.hpp" namespace "triqs::gf" : 

    cdef void test_gf_c "test_gf" ( gf_im_freq_c &)  
    cdef void test_block_c "test_block" ( gf_block_im_freq_c&)  
    cdef void test2_c "test2" ( matrix_view[double,COrder] &) except + 
    cdef void test3_c "test3" ( vector[matrix_view[double,COrder] ] &) except + 
    #cdef void test3_c "test3" ( vector[double] &)  

import numpy

def test_block (G) : 
    #cdef vector[gf_im_freq_c] v_c
    #for item in G:
    #    #v_c.push_back((<GFBloc_ImFreq_cython>item)._c)
    #    v_c.push_back(to_C_GFBloc_ImFreq(item))
    
    #cdef gf_block_view_c GG = make_gf_block_view_c (v_c)  
    test_block_c (as_gf_block_im_freq_c(G)) 
   
    #cdef GFBloc_ImFreq Res = 
    #Res._c = as_gf_block_im_freq_c(G)
    #return Res
     
    #test_gf_c (G[0]) 

# TEST OF AUTOMATIC CONVERSION
def test2():
    a = numpy.array([[1,2],[3,4]])
    #test2_c(a)

def test3():
    a = numpy.array([[1.0,2],[3,4]])
    w = [a,a]
    #test3_c(w)
    print a


