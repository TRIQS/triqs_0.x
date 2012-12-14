from dcomplex cimport * 
from arrays cimport *   
cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail "triqs::gf::local::tail_view"  :
        tail()
        tail(array_view[dcomplex,THREE,COrder] , int) except +
        matrix_view[dcomplex,COrder] operator()(int) except +
        array_view[dcomplex,THREE,COrder] data_view()
        int order_min()
        int order_max()
        size_t size() 
        size_t shape(int) 

    cdef tail operator +( tail &, tail &) except + 
    cdef tail operator -( tail &, tail &) except + 
    
    cdef tail operator *( tail&, tail &) 

    cdef tail operator *( double, tail &) except + 
    cdef tail operator *( tail &, double) except + 
    cdef tail operator /( double, tail &) except + 
    cdef tail operator /( tail &, double) except + 

    cdef tail operator *( dcomplex, tail &) except + 
    cdef tail operator *( tail &, dcomplex) except + 
    cdef tail operator /( dcomplex, tail &) except + 
    cdef tail operator /( tail &, dcomplex) except + 
    cdef tail inverse_c "inverse" ( tail &) except +

    cdef tail operator *( matrix_view[dcomplex,COrder] &, tail &) except + 
    cdef tail operator *( tail &, matrix_view[dcomplex,COrder]&) except + 

    cdef void h5_write (h5_group_or_file, char *, tail &)

