cdef extern from "triqs/gf/imtime.hpp" namespace "triqs::gf" : 
  
    cdef cppclass imtime_domain :
        double beta
        statistic_enum statistic
        imtime_domain ()

    cdef cppclass mesh_imtime "triqs::gf::linear_mesh<triqs::gf::imtime::domain_t>"  :
        mesh_imtime ()
        mesh_imtime (mesh_imtime &)
        imtime_domain & domain()
        long size()
        bint operator ==( mesh_imtime &)

    cdef mesh_imtime make_mesh_imtime "triqs::gf::imtime::make_mesh" (double beta, statistic_enum S, size_t Nmax)
    
    cdef cppclass gf_imtime "triqs::gf::gf_view<triqs::gf::imtime>" :
        gf_imtime()
        gf_imtime(gf_imtime &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imtime(mesh_imtime, array_view[double, THREE,COrder], tail, nothing, indices_2_t) #except +
        mesh_imtime mesh() 
        array_view[double, THREE,COrder] data_view()
        tail singularity_view() 
        indices_2_t indices()

    cdef gf_imtime operator +( gf_imtime &, gf_imtime &) except + 
    cdef gf_imtime operator -( gf_imtime &, gf_imtime &) except + 
    cdef gf_imtime operator *( gf_imtime &, gf_imtime &) except + 
    
    cdef gf_imtime operator *( double, gf_imtime &) except + 
    cdef gf_imtime operator *( gf_imtime &, double) except + 
    cdef gf_imtime operator /( gf_imtime &, double) except + 

    cdef gf_imtime operator *( matrix_view[double,COrder] &, gf_imtime &) except + 
    cdef gf_imtime operator *( gf_imtime &, matrix_view[double,COrder]&) except + 

cdef extern from "triqs/gf/imtime.hpp"  :

    cdef void h5_write (h5_group_or_file, char *, gf_imfreq &)
    cdef gf_imtime inverse_c "inverse"   (gf_imtime &)
    #cdef gf_imtime transpose_c (gf_imtime &)
    #cdef gf_imtime conjugate_c (gf_imtime &)
    
    cdef gf_imtime make_gf_imtime "triqs::gf::imtime::make_gf" (mesh_imtime, array_view[double, THREE,COrder], tail) except +
    cdef gf_imtime clone_gf_imtime "triqs::make_clone" (gf_imtime &) 

cdef extern from "triqs/utility/serialization.hpp"  :
    cdef std_string boost_serialize "triqs::serialize" (gf_imtime &) 
    cdef void boost_unserialize_into "triqs::deserialize_into_view" (std_string, gf_imtime &) 

# Python -> C
cdef gf_imtime  as_gf_imtime (g) except +  

# C -> Python 
cdef make_GfImTime ( gf_imtime x) 

###############  Blocks of Im Time #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_imtime "triqs::gf::gf_view<triqs::gf::block<triqs::gf::imtime> >" :
        gf_block_imtime()
        gf_imtime & operator [](int)
        discrete_mesh & mesh()

    cdef gf_block_imtime  make_gf_block_imtime "triqs::gf::block<triqs::gf::imtime>::make_gf_view" (  vector[gf_imtime] &) 

cdef gf_block_imtime  as_gf_block_imtime (G) except +
cdef make_BlockGfImTime (gf_block_imtime G) 


