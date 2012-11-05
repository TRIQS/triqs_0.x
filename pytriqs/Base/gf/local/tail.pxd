from dcomplex cimport * 
from arrays cimport *   
cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail_view_c "triqs::gf::local::tail_view"  :
        tail_view_c()
        tail_view_c(array_view[dcomplex,THREE,COrder] , int) #except +
        void rebind (tail_view_c &)

cdef class TailGF_c:
    cdef tail_view_c _c
 
