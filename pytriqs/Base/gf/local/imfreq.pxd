cdef extern from "triqs/gf/imfreq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass imfreq_domain :
        double beta
        statistic_enum statistic
        imfreq_domain ()

    cdef cppclass mesh_imfreq "triqs::gf::linear_mesh<triqs::gf::imfreq::domain_t>"  :
        mesh_imfreq ()
        mesh_imfreq (mesh_imfreq &)
        imfreq_domain & domain()
        long size()
        bint operator ==( mesh_imfreq &)

    cdef mesh_imfreq make_mesh_imfreq "triqs::gf::imfreq::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_imfreq "triqs::gf::gf_view<triqs::gf::imfreq>" :
        gf_imfreq()
        gf_imfreq(gf_imfreq &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imfreq(mesh_imfreq, array_view[dcomplex, THREE,COrder], tail, nothing) #except +
        mesh_imfreq mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail singularity_view() 

    cdef gf_imfreq operator +( gf_imfreq &, gf_imfreq &) except + 
    cdef gf_imfreq operator -( gf_imfreq &, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, gf_imfreq &) except + 
    
    cdef gf_imfreq operator *( dcomplex, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, dcomplex) except + 
    cdef gf_imfreq operator /( gf_imfreq &, dcomplex) except + 

    cdef gf_imfreq operator *( matrix_view[dcomplex,COrder] &, gf_imfreq &) except + 
    cdef gf_imfreq operator *( gf_imfreq &, matrix_view[dcomplex,COrder]&) except + 

cdef extern from "triqs/gf/imfreq.hpp"  :

    cdef gf_imfreq inverse   (gf_imfreq &)
    cdef gf_imfreq transpose (gf_imfreq &)
    cdef gf_imfreq conjugate (gf_imfreq &)
    
    cdef gf_imfreq make_gf_imfreq "triqs::gf::imfreq::make_gf" (mesh_imfreq, array_view[dcomplex, THREE,COrder], tail) except +
    cdef gf_imfreq clone_gf_imfreq "triqs::make_clone" (gf_imfreq &) 

cdef class _ImplGfLocal : 
    cdef object _myIndicesGFBlocL, _myIndicesGFBlocR, _Name, dtype, _IndicesR, _IndicesL

cdef class GfImFreq(_ImplGfLocal):
    cdef gf_imfreq _c

# Python -> C
cdef inline gf_imfreq  as_gf_imfreq (GfImFreq g) except +: 
    return g._c

# C -> Python 
cdef inline make_GfImFreq ( gf_imfreq x) : 
    #return GfImFreq(mesh = make_MeshImFreq(x.mesh()), data = x.data_view().to_python(), tail = make_TailGf(x.singularity_view()))
    return GfImFreq(C_Object = encapsulate (&x))

###############  Blocks of Im Freq #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_imfreq "triqs::gf::gf_view<triqs::gf::block<triqs::gf::imfreq> >" :
        gf_block_imfreq()
    
    cdef gf_block_imfreq  make_gf_block_imfreq "triqs::gf::block<triqs::gf::imfreq>::make_gf_view" (  vector[gf_imfreq] &) 

cdef inline gf_block_imfreq  as_gf_block_imfreq (G) except +:
    cdef vector[gf_imfreq] v_c
    for item in G:
        v_c.push_back((<GfImFreq>item)._c)
        v_c.push_back(as_gf_imfreq(item))
    return make_gf_block_imfreq (v_c) 


