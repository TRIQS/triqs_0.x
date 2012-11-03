#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
#from libcpp.vector cimport vector
cimport cython  
from view_proxy cimport *   
from gf_matsubara cimport *

cdef class TailGF_c:
    cdef view_proxy[tail_view] _c
    def __init__(self, a, int omin):
        self._c.rebind(  tail_view( array_view[dcomplex,THREE,COrder](a), omin))

include "./matsubara_freq.pyx"
include "./matsubara_time.pyx"  
