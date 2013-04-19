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

    cdef mesh_imfreq make_mesh_imfreq "triqs::gf::gf_factories<triqs::gf::imfreq>::make_mesh" (double beta, statistic_enum S, size_t n_max)
    
    cdef cppclass gf_imfreq "triqs::python_tools::cython_proxy<triqs::gf::gf_view<triqs::gf::imfreq> >" :
        gf_imfreq()
        gf_imfreq(gf_imfreq &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_imfreq(mesh_imfreq, array_view[dcomplex, THREE], tail, nothing) #except +
        void operator << (gf_imfreq &)
        mesh_imfreq mesh() 
        array_view[dcomplex, THREE] data_view()
        tail singularity_view() 

cdef extern from "triqs/gf/imfreq.hpp"  :
    cdef void h5_write (h5_group, char *, gf_imfreq &)

cdef extern from "triqs/utility/serialization.hpp"  :
    cdef std_string boost_serialize "triqs::serialize" (gf_imfreq &) 
    cdef void boost_unserialize_into "triqs::deserialize_into_view" (std_string, gf_imfreq &) 

# Python -> C
cdef gf_imfreq  as_gf_imfreq (g) except +  

# C -> Python 
cdef make_GfImFreq (gf_imfreq x, indices_pack=*, name=*)

###############  Blocks of Im Freq #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_imfreq "triqs::python_tools::cython_proxy<triqs::gf::gf_view<triqs::gf::block<triqs::gf::imfreq>>>" :
        gf_block_imfreq()
        gf_imfreq & operator [](int)
        discrete_mesh & mesh()

    cdef gf_block_imfreq  make_gf_block_imfreq "triqs::gf::make_gf_view<triqs::gf::block<triqs::gf::imfreq>>" (  vector[gf_imfreq] &) 

cdef gf_block_imfreq  as_gf_block_imfreq (G) except +
cdef make_BlockGfImFreq (gf_block_imfreq G, block_indices_pack=*, name=*)

