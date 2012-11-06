#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
cimport cython  

###############  Fourier  #########################

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_im_freq_c lazy_fourier          (gf_im_time_c & )
    gf_im_time_c lazy_inverse_fourier  (gf_im_freq_c & )

###############  Tail  #########################

cdef class TailGF_c:
    #cdef tail_view_c _c
    def __init__(self, a, int omin):
        self._c =  tail_view_c( array_view[dcomplex,THREE,COrder](a), omin)

###############  Frequences  #########################

# ----------- Domain  --------------------------
#cdef class DomainMatsubaraFrequency:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshImFreq: 
    #cdef matsubara_freq_mesh_c  _c

    def __init__(self, Beta, stat, int Nmax): 
        self._c =  matsubara_freq_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[matsubara_freq_mesh_c ] g = mesh_pt_generator[matsubara_freq_mesh_c ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- GF  --------------------------

cdef class GfImFreq_cython:
    #cdef object _mesh
    #cdef gf_im_freq_c _c
    
    def __init__(self, MeshImFreq m, data, TailGF_c tail):
        self._c =  gf_im_freq_c ( m._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing())
        self._mesh = m

    #def test(self) :
        #return toPy_gf_im_freq(self._c)
 
    def setFromFourierOf(self, gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( (<GfImTime_cython >gt)._c )

# ----------- Domain  --------------------------
#cdef class DomainMatsubaraTime:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshImTime: 
    #cdef matsubara_time_mesh_c _c
    
    def __init__(self, Beta, stat, int Nmax): 
        self._c = matsubara_time_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) :
        cdef mesh_pt_generator[matsubara_time_mesh_c] g = mesh_pt_generator[matsubara_time_mesh_c](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- The GF  --------------------------

cdef class GfImTime_cython:
    #cdef gf_im_time_c _c
    #cdef object _mesh
    def __init__(self, MeshImTime m, data, TailGF_c tail):
        self._c =  gf_im_time_c ( m._c, array_view[dcomplex,THREE,COrder](data), tail._c, nothing() )
        self._mesh = m
    
    def setFromInverseFourierOf(self,  gw) : 
        """Fills self with the Inverse Fourier transform of gw"""        
        self._c = lazy_inverse_fourier( (<GfImFreq_cython>gw)._c)

#---------------------- TEST FOR DEBUG ONLY ...------------------------------------------

cdef extern from "pytriqs/Base/gf/local/a.hpp" namespace "triqs::gf" : 

    cdef void test_gf_c "test_gf" ( gf_im_freq_c &)  
    cdef void test_block_c "test_block" ( gf_block_im_freq_c&)  
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
    test_block_c (as_gf_block_im_freq_c(G)) 
   
    #cdef GfImFreq Res = 
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


