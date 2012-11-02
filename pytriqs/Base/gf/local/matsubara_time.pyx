# ----------- C++  --------------------------
cdef extern from "triqs/gf/matsubara_time.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_time_domain :
        double beta
        statistic_enum statistic
        matsubara_time_domain ()

    cdef cppclass matsubara_time_mesh "triqs::gf::linear_mesh<triqs_gf_matsubara_time_domain_t>"  :
        double beta
        statistic_enum statistic
        matsubara_time_mesh ()
        matsubara_time_mesh (matsubara_time_mesh&)
        matsubara_time_domain & domain()
        long size()
    
    cdef matsubara_time_mesh matsubara_time_make_mesh "triqs::gf::matsubara_time::make_mesh" (double beta, statistic_enum S, size_t Nmax)

    cdef cppclass gf_view_time "triqs::gf::gf_view<triqs_gf_matsubara_time_desc>" :
        gf_view_time()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_view_time(matsubara_time_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view_time&)
        matsubara_time_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_view_time matsubara_time_make_gf "triqs::gf::matsubara_time::make_gf" (matsubara_time_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_view_time lazy_inverse_fourier (gf_view_freq & )

# ----------- Domain  --------------------------
#cdef class DomainMatsubaraTime:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshMatsubaraTime: 
    cdef matsubara_time_mesh _c
    
    def __init__(self, Beta, stat, int Nmax): 
        self._c = matsubara_time_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) :
        cdef mesh_pt_generator[matsubara_time_mesh] g = mesh_pt_generator[matsubara_time_mesh](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- The GF  --------------------------

cdef class GFBloc_ImTime_cython:
    cdef gf_view_time _c
    cdef object _mesh
    def __init__(self, MeshMatsubaraTime mesh, data, TailGF_c tail):
        self._c.rebind( gf_view_time ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c, nothing() ))
        self._mesh = mesh
    
    def setFromInverseFourierOf(self, GFBloc_ImFreq_cython gw) : 
        """Fills self with the Inverse Fourier transform of gw"""        
        self._c = lazy_inverse_fourier( gw._c)

