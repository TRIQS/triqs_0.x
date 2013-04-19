cdef extern from "triqs/gf/refreq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass refreq_domain :
        refreq_domain()

    cdef cppclass mesh_refreq "triqs::gf::linear_mesh<triqs::gf::refreq::domain_t>"  :
        mesh_refreq ()
        mesh_refreq (mesh_refreq &)
        refreq_domain & domain()
        double x_min()
        double x_max()
        long size()
        double kind()
        bint operator ==( mesh_refreq &)

    cdef mesh_refreq make_mesh_refreq "triqs::gf::gf_factories<triqs::gf::refreq>::make_mesh" (double omega_min, double omega_max, size_t n_freq, mesh_enum mk)

    cdef cppclass gf_refreq "triqs::python_tools::cython_proxy<triqs::gf::gf_view<triqs::gf::refreq>>" :
        gf_refreq()
        gf_refreq(gf_refreq &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_refreq(mesh_refreq, array_view[dcomplex, THREE], tail, nothing) #except +
        void operator << (gf_refreq &)
        mesh_refreq mesh() 
        array_view[dcomplex, THREE] data_view()
        tail singularity_view() 

cdef extern from "triqs/gf/refreq.hpp"  :
    cdef void h5_write (h5_group, char *, gf_refreq &)

cdef extern from "triqs/utility/serialization.hpp"  :
    cdef std_string boost_serialize "triqs::serialize" (gf_refreq &) 
    cdef void boost_unserialize_into "triqs::deserialize_into_view" (std_string, gf_refreq &) 

# Python -> C
cdef gf_refreq as_gf_refreq (g) except +

# C -> Python 
cdef make_GfReFreq (gf_refreq x, indices_pack=*, name=*)

###############  Blocks of Im Time #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_refreq "triqs::python_tools::cython_proxy<triqs::gf::gf_view<triqs::gf::block<triqs::gf::refreq>>>" :
        gf_block_refreq()
        gf_refreq & operator [](int)
        discrete_mesh & mesh()

    cdef gf_block_refreq  make_gf_block_refreq "triqs::gf::make_gf_view<triqs::gf::block<triqs::gf::refreq>>" (  vector[gf_refreq] &) 

cdef gf_block_refreq  as_gf_block_refreq (G) except +
cdef make_BlockGfReFreq (gf_block_refreq G, block_indices_pack=*, name=*)

