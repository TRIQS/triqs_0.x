from dcomplex cimport * 
from arrays cimport *   
from libcpp.vector cimport vector
from shared_ptr cimport *
from extractor cimport *

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
    cdef cppclass tail_view "triqs::gf::local::tail_view"  :
        tail_view()
        tail_view(array_view[dcomplex,THREE,COrder] , int) except +
        void rebind (tail_view &)
        tail_view operator +( tail_view &) 
        tail_view operator -( tail_view &)
        void print_me()
        array_view[dcomplex,THREE,COrder] data_view()
    cdef tail_view operator *( double, tail_view &) except + 
    cdef tail_view operator *( tail_view &, double) except + 
    cdef tail_view operator /( double, tail_view &) except + 
    cdef tail_view operator /( tail_view &, double) except + 

cdef class TailGF:
    cdef tail_view _c
    #cdef array_view[dcomplex,THREE,COrder] _ac

###############  IM FREQ #########################

cdef extern from "triqs/gf/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain :
        double beta
        statistic_enum statistic
        matsubara_freq_domain ()

    cdef cppclass matsubara_freq_mesh "triqs::gf::linear_mesh<triqs::gf::matsubara_freq::domain_t>"  :
        matsubara_freq_mesh ()
        matsubara_freq_mesh (matsubara_freq_mesh &)
        matsubara_freq_mesh(object,int)
        matsubara_freq_domain & domain()
        long size()
        bint operator ==( matsubara_freq_mesh &)

    cdef matsubara_freq_mesh matsubara_freq_make_mesh "triqs::gf::matsubara_freq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_imfreq "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_imfreq()
        gf_imfreq(gf_imfreq &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imfreq(matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        matsubara_freq_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_imfreq make_gf_imfreq "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

cdef class MeshImFreq: 
    cdef matsubara_freq_mesh  _c

cdef class GfImFreq_cython:
    cdef object _mesh
    cdef gf_imfreq _c

# Python -> C
cdef inline gf_imfreq  as_gf_imfreq (GfImFreq_cython g) except +: 
    return g._c

# C -> Python 
cdef inline make_GfImFreq ( gf_imfreq x) : 
    from gf_imfreq import GfImFreq 
    return GfImFreq(C_Object = encapsulate (&x))

###############  Blocks of Im Freq #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_imfreq "triqs::gf::gf_view<triqs::gf::block<triqs::gf::matsubara_freq> >" :
        gf_block_imfreq()
    
    cdef gf_block_imfreq  make_gf_block_imfreq "triqs::gf::block<triqs::gf::matsubara_freq>::make_gf_view" (  vector[gf_imfreq] &) 

cdef inline gf_block_imfreq  as_gf_block_imfreq (G) except +:
    cdef vector[gf_imfreq] v_c
    for item in G:
        v_c.push_back((<GfImFreq_cython>item)._c)
        v_c.push_back(as_gf_imfreq(item))
    return make_gf_block_imfreq (v_c) 



###############  IM TIME #########################

cdef extern from "triqs/gf/matsubara_time.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_time_domain :
        double beta
        statistic_enum statistic
        matsubara_time_domain ()

    cdef cppclass matsubara_time_mesh "triqs::gf::linear_mesh<triqs::gf::matsubara_time::domain_t>"  :
        matsubara_time_mesh ()
        matsubara_time_mesh (matsubara_time_mesh&)
        matsubara_time_domain & domain()
        long size()
    
    cdef matsubara_time_mesh matsubara_time_make_mesh "triqs::gf::matsubara_time::make_mesh" (double beta, statistic_enum S, size_t Nmax)

    cdef cppclass gf_imtime "triqs::gf::gf_view<triqs::gf::matsubara_time>" :
        gf_imtime()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imtime(matsubara_time_mesh, array_view[double, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_imtime&)
        matsubara_time_mesh mesh() 
        array_view[double, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_imtime matsubara_time_make_gf "triqs::gf::matsubara_time::make_gf" (matsubara_time_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_imtime "triqs::gf::gf_view<triqs::gf::block<triqs::gf::matsubara_time> >" :
        gf_block_imtime()
    
    cdef gf_block_imtime  make_gf_block_imtime "triqs::gf::block<triqs::gf::matsubara_time>::make_gf_view" (  vector[gf_imtime] &) 

cdef class MeshImTime:
    cdef matsubara_time_mesh _c

cdef class GfImTime_cython:
    cdef gf_imtime _c
    cdef object _mesh

cdef inline gf_imtime  as_gf_imtime (GfImTime_cython g) : 
    return g._c

cdef inline gf_block_imtime  as_gf_block_imtime (G):
    cdef vector[gf_imtime] v
    for item in G:
        v.push_back(as_gf_imtime(item))
    return make_gf_block_imtime (v) 


