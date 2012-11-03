from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
from gf_basic_tools cimport *
from tail cimport *

################   FREQ ################################

cdef extern from "triqs/gf/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain :
        double beta
        statistic_enum statistic
        matsubara_freq_domain ()

    cdef cppclass matsubara_freq_mesh "triqs::gf::linear_mesh<triqs::gf::matsubara_freq::domain_t>"  :
        double beta
        statistic_enum statistic
        matsubara_freq_mesh ()
        matsubara_freq_mesh (matsubara_freq_mesh &)
        matsubara_freq_domain & domain()
        long size()

    cdef matsubara_freq_mesh matsubara_freq_make_mesh "triqs::gf::matsubara_freq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_view_freq "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_view_freq()
        gf_view_freq(matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view_freq&)
        matsubara_freq_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_view_freq matsubara_freq_make_gf "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

################ TIME  ################################

cdef extern from "triqs/gf/matsubara_time.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_time_domain :
        double beta
        statistic_enum statistic
        matsubara_time_domain ()

    cdef cppclass matsubara_time_mesh "triqs::gf::linear_mesh<triqs::gf::matsubara_time::domain_t>"  :
        double beta
        statistic_enum statistic
        matsubara_time_mesh ()
        matsubara_time_mesh (matsubara_time_mesh&)
        matsubara_time_domain & domain()
        long size()
    
    cdef matsubara_time_mesh matsubara_time_make_mesh "triqs::gf::matsubara_time::make_mesh" (double beta, statistic_enum S, size_t Nmax)

    cdef cppclass gf_view_time "triqs::gf::gf_view<triqs::gf::matsubara_time>" :
        gf_view_time()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_view_time(matsubara_time_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view_time&)
        matsubara_time_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_view_time matsubara_time_make_gf "triqs::gf::matsubara_time::make_gf" (matsubara_time_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +


cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_view_freq lazy_fourier (gf_view_time & )
    gf_view_time lazy_inverse_fourier (gf_view_freq & )
