from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
#from libcpp.vector cimport vector
cimport cython
from arrays cimport *

cdef extern from "<complex>" namespace "std":
    cdef cppclass dcomplex "std::complex<double>":
         complex(double,double)

cdef extern from "triqs/gf/tools.hpp" namespace "triqs::gf" : 
    cdef enum statistic_enum "triqs::gf::statistic_enum" :
        Boson,Fermion

    cdef cppclass freq_infty: 
        freq_infty ()

    cdef cppclass mesh_pt_generator "triqs::gf::mesh_pt_generator" [MeshType] : 
        mesh_pt_generator( MeshType &, bool)
        complex to_point()
        mesh_pt_generator operator++()
        bint at_end()

cdef extern from "triqs/gf/descriptors/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain :
        double beta
        statistic_enum statistic
        matsubara_freq_domain ()

    cdef cppclass matsubara_freq_mesh " triqs::gf::matsubara_freq::mesh_t"  :
        double beta
        statistic_enum statistic
        matsubara_freq_mesh (double Beta, statistic_enum s, int Nmax) 
        matsubara_freq_domain & domain()
        long size()

cdef class DomainMatsubaraFrequency:
    pass

cdef class MeshMatsubaraFrequency_generator : 
    cdef mesh_pt_generator[matsubara_freq_mesh] * g
    def __init__(self,MeshMatsubaraFrequency m) : 
        self.g = new mesh_pt_generator[matsubara_freq_mesh](deref(m._c), True)
    def __next__(self) : 
        inc(deref(self.g))
        if self.g.at_end() : raise StopIteration
        return self.g.to_point()
    def __dealloc__(self): del self.g

cdef class MeshMatsubaraFrequency: 
    cdef matsubara_freq_mesh * _c
    def __init__(self, Beta, stat, int Nmax): self._c = new matsubara_freq_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    def __dealloc__(self): del self._c
    def __len__ (self) : return self._c.size()
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    def __iter__(self) :
        return MeshMatsubaraFrequency_generator(self)
    def __iter2__(self) :
        from math import pi
        for n in range(len(self)): 
            yield 1j*(2*n+1)*pi/self.Beta


cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail_view "triqs::gf::local::tail_view"  : 
        tail_view(array_view[dcomplex,THREE] , int) except +

cdef class TailGF_c:
    cdef tail_view * _c
    def __init__(self, a, int omin):
        self._c = new tail_view( array_view[dcomplex,THREE](a), omin)
    def __dealloc__(self):
        del self._c

cdef extern from "triqs/gf/gf.hpp" namespace "triqs::gf" : 
    cdef cppclass gf_view "triqs::gf::gf_view<triqs::gf::matsubara_freq>"  : 
        gf_view(matsubara_freq_mesh, array_view[dcomplex, THREE], tail_view) except +

cdef class GFBloc_ImFreq_c:
    cdef gf_view *_c
    def __init__(self, MeshMatsubaraFrequency mesh, data, TailGF_c tail):
        self._c = new gf_view ( deref(mesh._c), array_view[dcomplex,THREE](data), deref(tail._c) )
    def __dealloc__(self):
        del self._c
