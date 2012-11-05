from arrays cimport *   
from shared_ptr cimport *
from dcomplex cimport * 
from gf_basic_tools cimport *
from tail cimport *

################ GAL  ################################
# DOES NOT WORK YET ...
# cython has a pb with ---> get gf::DESC><domain_t> etc...
# moreover, the domain and meshes will have non generic parts...
# Is it a good idea to write something general ???
cdef extern from "triqs/gf/gf.hpp"  : 
    cdef cppclass gf_domain "triqs::gf::DESC::domain_t" [DESC] : 
        gf_domain()
        # NOT always present : error will detected at C++ compilation
        double beta
        statistic_enum statistic
 
    cdef cppclass linear_mesh "triqs::gf::linear_mesh<triqs::gf::DESC::domain_t>" [DESC] :
        linear_mesh ()
        linear_mesh (linear_mesh[DESC]&)
        gf_domain[DESC] & domain()
        long size()

    cdef cppclass gf_mesh "triqs::gf::DESC::mesh_t" [DESC] : 
        gf_mesh()
        gf_mesh (gf_mesh[DESC]&)
        gf_domain[DESC] & domain()
        long size()

    cdef cppclass gf_view "triqs::gf::gf_view<triqs::gf::DESC>" [DESC]:
        gf_view()
        gf_view(gf_mesh[DESC], array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view&)
        gf_mesh[DESC] mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

