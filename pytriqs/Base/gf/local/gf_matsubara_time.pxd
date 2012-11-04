from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
from gf_basic_tools cimport *
from tail cimport *

cdef extern from "triqs/gf/matsubara_time.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_time_domain_c :
        double beta
        statistic_enum statistic
        matsubara_time_domain_c ()

    cdef cppclass matsubara_time_mesh_c "triqs::gf::linear_mesh<triqs::gf::matsubara_time::domain_t>"  :
        matsubara_time_mesh_c ()
        matsubara_time_mesh_c (matsubara_time_mesh_c&)
        matsubara_time_domain_c & domain()
        long size()
    
    cdef matsubara_time_mesh_c matsubara_time_make_mesh "triqs::gf::matsubara_time::make_mesh" (double beta, statistic_enum S, size_t Nmax)

    cdef cppclass gf_view_time_c "triqs::gf::gf_view<triqs::gf::matsubara_time>" :
        gf_view_time_c()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_view_time_c(matsubara_time_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c, nothing) #except +
        void rebind( gf_view_time_c&)
        matsubara_time_mesh_c mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view_c auxiliary_data_view() 

    cdef gf_view_time_c matsubara_time_make_gf "triqs::gf::matsubara_time::make_gf" (matsubara_time_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c) except +

