from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
from gf_basic_tools cimport *
from tail cimport *

cdef extern from "triqs/gf/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain :
        double beta
        statistic_enum statistic
        matsubara_freq_domain ()

    cdef cppclass matsubara_freq_mesh "triqs::gf::linear_mesh<triqs::gf::matsubara_freq::domain_t>"  :
        matsubara_freq_mesh ()
        matsubara_freq_mesh (matsubara_freq_mesh &)
        matsubara_freq_domain & domain()
        long size()
   
    cdef matsubara_freq_mesh matsubara_freq_make_mesh "triqs::gf::matsubara_freq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_view_freq "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_view_freq()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_view_freq(matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view_freq&)
        matsubara_freq_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_view_freq matsubara_freq_make_gf "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

