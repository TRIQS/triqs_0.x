from dcomplex cimport * 
from arrays cimport *   
cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail "triqs::python_tools::cython_proxy<triqs::gf::local::tail_view>"  :
        tail()
        tail(array_view[dcomplex,THREE], int, array_view[long,TWO]) except +
        matrix_view[dcomplex] operator()(int) except +
        array_view[dcomplex,THREE] data_view()
        array_view[long,TWO] mask_view()
        void operator << (tail &)
        long order_min()
        long order_max()
        size_t size() 
        size_t shape(int) 

    cdef tail operator +( tail &, tail &) except + 
    cdef tail operator -( tail &, tail &) except + 
    
    cdef tail operator *( tail&, tail &) except +

    cdef tail operator *( double, tail &) except + 
    cdef tail operator *( tail &, double) except + 
    cdef tail operator /( double, tail &) except + 
    cdef tail operator /( tail &, double) except + 

    cdef tail operator *( dcomplex, tail &) except + 
    cdef tail operator *( tail &, dcomplex) except + 
    cdef tail operator /( dcomplex, tail &) except + 
    cdef tail operator /( tail &, dcomplex) except + 
    cdef tail inverse_c "inverse" ( tail &) except +

    cdef tail operator *( matrix_view[dcomplex] &, tail &) except + 
    cdef tail operator *( tail &, matrix_view[dcomplex]&) except + 

    cdef void h5_write (h5_group_or_file, char *, tail &)

cdef extern from "triqs/utility/serialization.hpp"  :
    cdef std_string boost_serialize "triqs::serialize" (tail &)
    cdef void boost_unserialize_into "triqs::deserialize_into_view" (std_string, tail &)
