from dcomplex cimport * 
from arrays cimport *   
cdef extern from "triqs/gf/local/tail.hpp" : 
    cdef cppclass tail_view "triqs::gf::local::tail_view"  :
        tail_view()
        tail_view(array_view[dcomplex,THREE,COrder] , int) #except +
        void rebind (tail_view &)


