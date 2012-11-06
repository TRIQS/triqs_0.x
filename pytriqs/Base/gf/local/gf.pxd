from dcomplex cimport * 
from arrays cimport *   
from libcpp.vector cimport vector
from shared_ptr cimport *

# -------------------- Some generic tools -------------------------------
 
cdef extern from "triqs/gf/tools.hpp" namespace "triqs::gf" : 
    cdef enum statistic_enum "triqs::gf::statistic_enum" :
        Boson,Fermion 
    
    cdef cppclass nothing : 
        nothing ()

    cdef cppclass freq_infty: 
        freq_infty ()

    cdef cppclass mesh_pt_generator "triqs::gf::mesh_pt_generator" [MeshType] :  
        mesh_pt_generator( MeshType * )
        mesh_pt_generator()
        complex to_point()
        mesh_pt_generator operator++()
        bint at_end()
        void increment() 
 
###############  TAIL #########################

cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail_view_c "triqs::gf::local::tail_view"  :
        tail_view_c()
        tail_view_c(array_view[dcomplex,THREE,COrder] , int) #except +
        void rebind (tail_view_c &)

cdef class TailGF_c:
    cdef tail_view_c _c

###############  IM FREQ #########################

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
    
    cdef cppclass gf_im_freq_c "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_im_freq_c()
        gf_im_freq_c(gf_im_freq_c &)
        gf_im_freq_c(object,int)
        object address_as_opaque_python_object()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_im_freq_c(matsubara_freq_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c, nothing) #except +
        void rebind( gf_im_freq_c&)
        matsubara_freq_mesh_c mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view_c auxiliary_data_view() 

    cdef gf_im_freq_c matsubara_freq_make_gf "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c) except +

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_im_freq_c "triqs::gf::gf_view<triqs::gf::block<matsubara_freq> >" :
        gf_block_im_freq_c()
    
    cdef gf_block_im_freq_c  make_gf_block_im_freq_c "triqs::gf::block<matsubara_freq>::make_gf_view" (  vector[gf_im_freq_c] &) 

cdef class MeshImFreq: 
    cdef matsubara_freq_mesh_c  _c

cdef class GfImFreq_cython:
    cdef object _mesh
    cdef gf_im_freq_c _c

cdef inline gf_im_freq_c  as_gf_im_freq_c (GfImFreq_cython g) except +: 
    return g._c

cdef inline gf_block_im_freq_c  as_gf_block_im_freq_c (G) except +:
    cdef vector[gf_im_freq_c] v_c
    for item in G:
        v_c.push_back((<GfImFreq_cython>item)._c)
        v_c.push_back(as_gf_im_freq_c(item))
    return make_gf_block_im_freq_c (v_c) 

from gf_im_freq import GfImFreq 

#cdef inline MeshImFreq_from_c ( matsubara_freq_mesh_c & c) :
#    return MeshImFreq(c.beta(),c.Nmax

#cdef inline GfImFreq_from_c  ( gf_im_freq_c& c) : 
#    return GfImFreq( MeshImFreq_from_c(c.mesh()), c.data_view(), TailGf_from_c (c.tail()), nothing())


###############  IM TIME #########################

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

    cdef cppclass gf_im_time_c "triqs::gf::gf_view<triqs::gf::matsubara_time>" :
        gf_im_time_c()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_im_time_c(matsubara_time_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c, nothing) #except +
        void rebind( gf_im_time_c&)
        matsubara_time_mesh_c mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view_c auxiliary_data_view() 

    cdef gf_im_time_c matsubara_time_make_gf "triqs::gf::matsubara_time::make_gf" (matsubara_time_mesh_c, array_view[dcomplex, THREE,COrder], tail_view_c) except +

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_im_time_c "triqs::gf::gf_view<triqs::gf::block<matsubara_time> >" :
        gf_block_im_time_c()
    
    cdef gf_block_im_time_c  make_gf_block_im_time_c "triqs::gf::block<matsubara_time>::make_gf_view" (  vector[gf_im_time_c] &) 

cdef class MeshImTime:
    cdef matsubara_time_mesh_c _c

cdef class GfImTime_cython:
    cdef gf_im_time_c _c
    cdef object _mesh

cdef inline gf_im_time_c  as_gf_im_time_c (GfImTime_cython g) : 
    return g._c

cdef inline gf_block_im_time_c  as_gf_block_im_time_c (G):
    cdef vector[gf_im_time_c] v_c
    for item in G:
        v_c.push_back(as_gf_im_time_c(item))
    return make_gf_block_im_time_c (v_c) 


