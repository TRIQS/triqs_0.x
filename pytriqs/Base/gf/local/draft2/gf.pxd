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

cdef class TailGf_cython:
    cdef tail_view _c
    #cdef array_view[dcomplex,THREE,COrder] _ac

# C -> Python 
cdef inline make_TailGf ( tail_view x) : 
    from gf_imfreq import GfImFreq 
    return None
    #return TailGf(mesh = make_MeshImFreq(x.mesh()), data = x.data_view().to_python(), tail = make_TailGf(x.singularity_view()))
    #return GfImFreq(C_Object = encapsulate (&x))


###############  IM FREQ #########################

cdef extern from "triqs/gf/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass imfreq_domain :
        double beta
        statistic_enum statistic
        imfreq_domain ()

    cdef cppclass mesh_imfreq "triqs::gf::linear_mesh<triqs::gf::matsubara_freq::domain_t>"  :
        mesh_imfreq ()
        mesh_imfreq (mesh_imfreq &)
        imfreq_domain & domain()
        long size()
        bint operator ==( mesh_imfreq &)

    cdef mesh_imfreq make_mesh_imfreq "triqs::gf::matsubara_freq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_imfreq "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_imfreq()
        gf_imfreq(gf_imfreq &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imfreq(mesh_imfreq, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        mesh_imfreq mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view singularity_view() 

    cdef gf_imfreq operator +( gf_imfreq &, gf_imfreq &) except + 
    cdef gf_imfreq operator -( gf_imfreq &, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, gf_imfreq &) except + 

    cdef gf_imfreq operator *( dcomplex, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, dcomplex) except + 
    cdef gf_imfreq operator /( gf_imfreq &, dcomplex) except + 

    cdef gf_imfreq operator *( matrix_view[dcomplex,COrder] &, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, matrix_view[dcomplex,COrder]&) except + 
    cdef gf_imfreq operator *( matrix[dcomplex,COrder] &, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, matrix[dcomplex,COrder]&) except + 

    cdef gf_imfreq inverse   (gf_imfreq &)
    cdef gf_imfreq transpose (gf_imfreq &)
    cdef gf_imfreq conjugate (gf_imfreq &)
    
    cdef gf_imfreq make_gf_imfreq "triqs::gf::matsubara_freq::make_gf" (mesh_imfreq, array_view[dcomplex, THREE,COrder], tail_view) except +
    cdef gf_imfreq clone_gf_imfreq "triqs::make_clone" (gf_imfreq &) 

###############  Cython classes #########################

cdef class MeshImFreq: 
    cdef mesh_imfreq  _c

# C -> Python 
cdef inline make_MeshImFreq ( mesh_imfreq x) :
    return MeshImFreq(C_Object = encapsulate (&x))

cdef class Gf_general : 
    pass

cdef class GfImFreq_cython (Gf_general):
    cdef gf_imfreq _c
    cdef MeshImFreq _mesh
    cdef TailGf_cython _tail

# Python -> C
cdef inline gf_imfreq  as_gf_imfreq (GfImFreq_cython g) except +: 
    return g._c

# C -> Python 
cdef inline make_GfImFreq ( gf_imfreq x) : 
    from gf_imfreq import GfImFreq 
    #return GfImFreq(mesh = make_MeshImFreq(x.mesh()), data = x.data_view().to_python(), tail = make_TailGf(x.singularity_view()))
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
  
    cdef cppclass imtime_domain :
        double beta
        statistic_enum statistic
        imtime_domain ()

    cdef cppclass mesh_imtime "triqs::gf::linear_mesh<triqs::gf::matsubara_time::domain_t>"  :
        mesh_imtime ()
        mesh_imtime (mesh_imtime&)
        imtime_domain & domain()
        long size()
    
    cdef mesh_imtime make_mesh_imtime "triqs::gf::matsubara_time::make_mesh" (double beta, statistic_enum S, size_t Nmax)

    cdef cppclass gf_imtime "triqs::gf::gf_view<triqs::gf::matsubara_time>" :
        gf_imtime()
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imtime(mesh_imtime, array_view[double, THREE,COrder], tail_view, nothing) #except +
        mesh_imtime mesh() 
        array_view[double, THREE,COrder] data_view()
        tail_view singularity_view() 

    cdef gf_imtime imtime_make_gf "triqs::gf::matsubara_time::make_gf" (mesh_imtime, array_view[dcomplex, THREE,COrder], tail_view) except +

cdef class MeshImTime:
    cdef mesh_imtime _c

cdef class GfImTime_cython:
    cdef gf_imtime _c
    cdef object _mesh

cdef inline gf_imtime  as_gf_imtime (GfImTime_cython g) : 
    return g._c

###############  Blocks of Im Time #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_imtime "triqs::gf::gf_view<triqs::gf::block<triqs::gf::matsubara_time> >" :
        gf_block_imtime()
    
    cdef gf_block_imtime  make_gf_block_imtime "triqs::gf::block<triqs::gf::matsubara_time>::make_gf_view" (  vector[gf_imtime] &) 

cdef inline gf_block_imtime  as_gf_block_imtime (G):
    cdef vector[gf_imtime] v
    for item in G:
        v.push_back(as_gf_imtime(item))
    return make_gf_block_imtime (v) 


