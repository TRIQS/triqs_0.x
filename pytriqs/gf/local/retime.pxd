cdef extern from "triqs/gf/retime.hpp" namespace "triqs::gf" : 
  
    cdef cppclass retime_domain :
        retime_domain()

    cdef cppclass mesh_retime "triqs::gf::linear_mesh<triqs::gf::retime::domain_t>"  :
        mesh_retime ()
        mesh_retime (mesh_retime &)
        retime_domain & domain()
        double x_min()
        double x_max()
        long size()
        double kind()
        bint operator ==( mesh_retime &)

    cdef mesh_retime make_mesh_retime "triqs::gf::gf_factories<triqs::gf::retime>::make_mesh" (double t_min, double t_max, size_t n_freq, mesh_enum mk)

    cdef cppclass gf_retime "triqs::python_tools::cython_proxy<triqs::gf::gf_view<triqs::gf::retime>>" :
        gf_retime()
        gf_retime(gf_retime &)
        # The constructor must be no_except, or the cython code won't be correct...
        gf_retime(mesh_retime, array_view[dcomplex, THREE], tail, nothing) #except +
        void operator << (gf_retime &)
        mesh_retime mesh() 
        array_view[dcomplex, THREE] data_view()
        tail singularity_view() 

cdef extern from "triqs/gf/retime.hpp"  :
    cdef void h5_write (h5_group, char *, gf_retime &)

cdef extern from "triqs/utility/serialization.hpp"  :
    cdef std_string boost_serialize "triqs::serialize" (gf_retime &) 
    cdef void boost_unserialize_into "triqs::deserialize_into_view" (std_string, gf_retime &) 

# Python -> C
cdef gf_retime as_gf_retime (g) except +

# C -> Python 
cdef make_GfReTime (gf_retime x, indices_pack=*, name=*)

###############  Blocks of Im Time #########################

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 

    cdef cppclass gf_block_retime "triqs::python_tools::cython_proxy<triqs::gf::gf_view<triqs::gf::block<triqs::gf::retime>>>" :
        gf_block_retime()
        gf_retime & operator [](int)
        discrete_mesh & mesh()

    cdef gf_block_retime  make_gf_block_retime "triqs::gf::make_gf_view<triqs::gf::block<triqs::gf::retime>>" (  vector[gf_retime] &) 

cdef gf_block_retime  as_gf_block_retime (G) except +
cdef make_BlockGfReTime (gf_block_retime G, block_indices_pack=*, name=*)

