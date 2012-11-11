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

cdef class TailGf_cython:
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

# ----------- Mesh  --------------------------
cdef class MeshImFreq: 
    #cdef mesh_imfreq  _c

    def __init__(self, Beta, stat, int Nmax): 
        self._c =  make_mesh_imfreq(Beta,{'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[mesh_imfreq ] g = mesh_pt_generator[mesh_imfreq ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

    def __richcmp__(MeshImFreq self, MeshImFreq other,int op) : 
        if op ==2 : # ==
            return self._c == other._c

# ----------- GF  --------------------------


cdef class Gf_general : 
    property tail : 
        def __get__(self): return self._tail
 
    def test_1(self) : 
        return self._c.data_view().to_python()


cdef class GfImFreq_cython(Gf_general):
    #cdef gf_imfreq _c
    #cdef MeshImFreq _mesh
    #cdef TailGf _tail
    def __init__(self, MeshImFreq mesh, data, TailGf_cython tail):
            self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing())
            self._mesh = mesh
            self._tail = tail

    # properties privatise element....
    property mesh : 
        """Mesh"""
        def __get__(self): return self._mesh
    
    property tail : 
        def __get__(self): return self._tail

    property _data_raw : 
        """Access to the data array"""
    
        def __get__(self) : 
            return self._c.data_view().to_python()

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            a[:,:,:] = value
    
    property _data : 
        """Access to the data array"""
    
        def __get__(self) : 
            return ArrayViewWithIndexConverter(self._c.data_view().to_python(), self._indL, self._indR, None)

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            a[:,:,:] = value
 
    def __div__(self, y) : 
        print "div " , self, y

    def __init2__(self, **d):

        c_obj = d.pop('C_Object', None)
        if c_obj :
            assert d == {}
            self._c = extractor [gf_imfreq] (c_obj) () 
            self._mesh = make_MeshImFreq (self._c.mesh())
        else : 
            def deleg (MeshImFreq mesh, data, TailGf_cython tail):
                self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing())
                self._mesh = mesh
            deleg(**d)

    def test(self) :
        return make_GfImFreq(self._c)
 
    def setFromFourierOf(self, gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( (<GfImTime_cython >gt)._c )

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (GfImFreq_cython)


def my_fnt( GfImFreq_cython g) : 
    print "ok", g._c.mesh().domain().beta

fd ={ 'g1' : my_fnt } 

include "_gf_imfreq.pyx"

# ----------- Mesh  --------------------------
cdef class MeshImTime: 
    #cdef mesh_imtime_c _c
    
    def __init__(self, Beta, stat, int Nmax): 
        self._c = make_mesh_imtime(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) :
        cdef mesh_pt_generator[mesh_imtime] g = mesh_pt_generator[mesh_imtime](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()
 
# ----------- The GF  --------------------------

cdef class GfImTime_cython:
    #cdef gf_imtime _c
    #cdef object _mesh
    def __init__(self, MeshImTime m, data, TailGf_cython tail):
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


