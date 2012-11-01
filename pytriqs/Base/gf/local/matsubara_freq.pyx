# ----------- C++  --------------------------
cdef extern from "triqs/gf/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain :
        double beta
        statistic_enum statistic
        matsubara_freq_domain ()

    cdef cppclass matsubara_freq_mesh "triqs::gf::linear_mesh<triqs_gf_matsubara_freq_domain_t>"  :
        double beta
        statistic_enum statistic
        matsubara_freq_mesh ()
        matsubara_freq_mesh (matsubara_freq_mesh &)
        matsubara_freq_domain & domain()
        long size()

    cdef matsubara_freq_mesh matsubara_freq_make_mesh "triqs::gf::matsubara_freq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    #cdef cppclass gf_view_freq "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
    cdef cppclass gf_view_freq "triqs::gf::gf_view<triqs_gf_matsubara_freq_desc>" :
        gf_view_freq()
        gf_view_freq(matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view_freq&)
        matsubara_freq_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_view_freq matsubara_freq_make_gf "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_view_freq lazy_fourier (gf_view_time & )

# ----------- Domain  --------------------------
#cdef class DomainMatsubaraFrequency:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshMatsubaraFrequency: 
    cdef matsubara_freq_mesh  _c
    
    def __init__(self, Beta, stat, int Nmax): 
        self._c =  matsubara_freq_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[matsubara_freq_mesh] g = mesh_pt_generator[matsubara_freq_mesh](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- GF  --------------------------

cdef class GFBloc_ImFreq_cython:
    cdef gf_view_freq _c
    cdef object _mesh
    
    def __init__(self, MeshMatsubaraFrequency mesh, data, TailGF_c tail):
            self._c.rebind( gf_view_freq( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing()))
            self._mesh = mesh
    
    def setFromFourierOf(self, GFBloc_ImTime_cython gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( gt._c)


