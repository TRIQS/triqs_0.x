#!python
#cython: embedsignature=True
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
cimport cython  
from pytriqs.Base.GF_Local.ArrayViewWithIndexConverter import ArrayViewWithIndexConverter

###############  Fourier  #########################

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_imfreq lazy_fourier          (gf_imtime & )
    gf_imtime lazy_inverse_fourier  (gf_imfreq & )

###############  Tail  #########################

cdef class TailGF:
    #cdef tail_view _c
    def __init__(self, a, int omin):
        self._c =  tail_view( array_view[dcomplex,THREE,COrder](a), omin)
    def invert(self) :
        self._c = 1.0/self._c
    
    def print_me(self) : 
        self._c.print_me()

    def __fill_a( self,a, value) : 
        if hasattr(value,'shape') : 
            if a.shape[:2] != value.shape[:2] : 
                raise RuntimeError, "shape mismatch"
            m = min(value.shape[2],a.shape[2]) 
            if a.shape[2] > value.shape[2] : 
                a[:,:,m:] =0
            a[:,:, 0:m] = value[:,:, 0:m]
        else : 
            a = value

    property _data_raw : 
        """Access to the data array"""
    
        def __get__(self) : 
            return self._c.data_view().to_python()

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            self.__fill_a(a,value)   
    
    property _data : 
        """Access to the data array"""
    
        def __get__(self) : 
            return ArrayViewWithIndexConverter(self._c.data_view().to_python(), self._indL, self._indR, None)

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            self.__fill_a(a,value)   
 
###############  Frequences  #########################

# ----------- Domain  --------------------------
#cdef class DomainMatsubaraFrequency:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshImFreq: 
    #cdef matsubara_freq_mesh  _c

    def __init__(self, Beta, stat, int Nmax): 
        self._c =  matsubara_freq_make_mesh(Beta,{'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[matsubara_freq_mesh ] g = mesh_pt_generator[matsubara_freq_mesh ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

    def __richcmp__(MeshImFreq self, MeshImFreq other,int op) : 
        if op ==2 : # ==
            return self._c == other._c

# C -> Python 
cdef inline make_MeshImFreq ( matsubara_freq_mesh x) :
    return MeshImFreq(C_Object = encapsulate (&x))

# ----------- GF  --------------------------

cdef class GfImFreq_cython:
    #cdef object _mesh
    #cdef gf_imfreq _c
    
    def __init__(self, **d):

        c_obj = d.pop('C_Object', None)
        if c_obj :
            assert d == {}
            self._c = extractor [gf_imfreq] (c_obj) () 
            self._mesh = make_MeshImFreq (self._c.mesh())
        else : 
            def deleg (MeshImFreq mesh, data, TailGF tail):
                self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing())
                self._mesh = mesh
            deleg(**d)

    def test(self) :
        return make_GfImFreq(self._c)
 
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
        cdef mesh_pt_generator[matsubara_time_mesh] g = mesh_pt_generator[matsubara_time_mesh](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()
 
# ----------- The GF  --------------------------

cdef class GfImTime_cython:
    #cdef gf_imtime _c
    #cdef object _mesh
    def __init__(self, MeshImTime m, data, TailGF tail):
        self._c =  gf_imtime ( m._c, array_view[double,THREE,COrder](data), tail._c, nothing() )
        self._mesh = m
    
    def setFromInverseFourierOf(self,  gw) : 
        """Fills self with the Inverse Fourier transform of gw"""        
        self._c = lazy_inverse_fourier( (<GfImFreq_cython>gw)._c)

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


