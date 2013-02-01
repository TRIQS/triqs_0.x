cdef extern from "triqs/gf/legendre.hpp" namespace "triqs::gf" : 
  
    cdef cppclass legendre_domain :
        size_t Nmax
        double beta
        statistic_enum statistic
        legendre_domain ()

    cdef cppclass mesh_legendre "triqs::gf::discrete_mesh<triqs::gf::legendre::domain_t>":
        mesh_legendre ()
        mesh_legendre (mesh_legendre &)
        legendre_domain & domain()
        long size()
        bint operator == (mesh_legendre &)

    cdef mesh_legendre make_mesh_legendre "triqs::gf::legendre::make_mesh" (double beta, statistic_enum S, size_t n_leg)

    cdef cppclass gf_legendre "triqs::gf::gf_view<triqs::gf::legendre>" :
        gf_legendre()
        gf_legendre(gf_legendre &)
        gf_legendre(mesh_legendre, array_view[dcomplex, THREE,COrder], nothing, nothing, indices_2_t) #except +
        mesh_legendre mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        indices_2_t indices()

cdef extern from "triqs/gf/legendre.hpp"  :

    cdef void h5_write (h5_group_or_file, char *, gf_legendre &)

cdef extern from "triqs/utility/serialization.hpp"  :
    cdef std_string boost_serialize "triqs::serialize" (gf_legendre &) 
    cdef void boost_unserialize_into "triqs::deserialize_into_view" (std_string, gf_legendre &) 

# Python -> C
cdef gf_legendre as_gf_legendre(g) except +  

# C -> Python 
cdef make_GfLegendre(gf_legendre x)

###############  Blocks of Im Time #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_legendre "triqs::gf::gf_view<triqs::gf::block<triqs::gf::legendre> >" :
        gf_block_legendre()
        gf_legendre & operator [](int)
        discrete_mesh & mesh()

    cdef gf_block_legendre  make_gf_block_legendre "triqs::gf::block<triqs::gf::legendre>::make_gf_view" (  vector[gf_legendre] &) 

cdef gf_block_legendre  as_gf_block_legendre (G) except +
cdef make_BlockGfLegendre (gf_block_legendre G) 


