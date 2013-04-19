from dcomplex cimport * 
from arrays cimport *   
from libcpp.vector cimport vector
from libcpp.string cimport string as std_string
from shared_ptr cimport *
from extractor cimport *
from h5 cimport *

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

    cdef cppclass indices_2_t : 
        indices_2_t( )
        indices_2_t( vector[vector[std_string]] &)
        vector[std_string] & operator[](int)
        vector[vector[std_string]] & operator()()
        bint same() 

cdef extern from "triqs/gf/meshes/linear.hpp" namespace "triqs::gf" :

    cdef enum mesh_enum "triqs::gf::mesh_kind":
        half_bins, full_bins, without_last

cdef extern from "triqs/gf/block.hpp" namespace "triqs::gf" : 
  
    cdef cppclass discrete_domain :
        discrete_domain ()
        discrete_domain (vector[std_string] &)
        vector[std_string] & names()

    cdef cppclass discrete_mesh :
        discrete_mesh ()
        discrete_mesh (discrete_mesh &)
        discrete_domain & domain()
        long size()
        bint operator ==( discrete_mesh &)

include "tail.pxd"
include "imfreq.pxd"
include "imtime.pxd"
include "refreq.pxd"
include "retime.pxd"
include "legendre.pxd"

