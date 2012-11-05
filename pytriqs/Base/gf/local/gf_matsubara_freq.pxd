from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
from gf_basic_tools cimport *
from tail cimport *

cdef extern from "triqs/gf/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain_c :
        double beta
        statistic_enum statistic
        matsubara_freq_domain_c ()

    cdef cppclass matsubara_freq_mesh_c "triqs::gf::linear_mesh<triqs::gf::matsubara_freq::domain_t>"  :
        matsubara_freq_mesh_c ()
        matsubara_freq_mesh_c (matsubara_freq_mesh_c &)
        matsubara_freq_domain_c & domain()
        long size()
   
    cdef matsubara_freq_mesh_c matsubara_freq_make_mesh "triqs::gf::matsubara_freq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_view_freq_c "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_view_freq_c()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_view_freq_c(matsubara_freq_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c, nothing) #except +
        void rebind( gf_view_freq_c&)
        matsubara_freq_mesh_c mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view_c auxiliary_data_view() 

    cdef gf_view_freq_c matsubara_freq_make_gf "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c) except +

cdef class MeshMatsubaraFrequency: 
    cdef matsubara_freq_mesh_c  _c

cdef class GFBloc_ImFreq_cython:
    cdef object _mesh
    cdef gf_view_freq_c _c

#cdef inline gf_view_freq_c   to_gf_view_freq_c (GFBloc_ImFreq_cython g) : 
#    return g._c

#cdef inline gf_block_view_c  to_gf_block_view_c (G):
#    cdef vector[gf_view_freq_c] v_c
#    for item in G:
#        v_c.push_back(to_gf_view_freq_c(item))
#    return make_gf_block_view_c (v_c) 



