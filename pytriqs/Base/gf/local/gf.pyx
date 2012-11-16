#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
cimport cython  
import numpy
import string
import warnings
from GF import GF
from math import pi
from h5 cimport *

include "h5.pxd"
include "fourier.pxd"
include "tail.pyx"
include "common.pyx"
include "imfreq.pyx"
include "imtime.pyx"

#---------------------- TEST FOR DEBUG ONLY ...------------------------------------------

cdef extern from "pytriqs/Base/gf/local/a.hpp" namespace "triqs::gf" : 

    cdef void test_gf_c "test_gf" ( gf_imfreq &)  
    cdef void test_block_c "test_block" ( gf_block_imfreq&)  
    cdef void test2_c "test2" ( matrix_view[double,COrder] &) except + 
    cdef void test3_c "test3" ( vector[matrix_view[double,COrder] ] &) except + 
    #cdef void test3_c "test3" ( vector[double] &)  

import numpy

def test_block (G) : 
    #cdef vector[gf_im_freq] v_c
    #for item in G:
    #    #v_c.push_back((<GfImFreq_cython>item)._c)
    #    v_c.push_back(to_C_GfImFreq(item))
    
    #cdef gf_block_view_c GG = make_gf_block_view_c (v_c)  
    test_block_c (as_gf_block_imfreq(G)) 
   
    #cdef GfImFreq Res = 
    #Res._c = as_gf_block_imfreq(G)
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


