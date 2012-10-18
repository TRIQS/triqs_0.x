cdef extern from "<triqs/arrays.hpp>" namespace "triqs::arrays" : 
    ctypedef int ONE   "1" 
    ctypedef int TWO   "2" 
    ctypedef int THREE "3" 
    ctypedef int FOUR  "4" 
    ctypedef int FIVE  "5" 
    ctypedef int SIX   "6" 

    cdef cppclass COrder "triqs::arrays::Option::C":
        COrder()
    cdef cppclass FortranOrder "triqs::arrays::Option::Fortran" :
        FortranOrder()

    # what about the fortran arrays ?
    cdef cppclass array_view "triqs::arrays::array_view" [T,R,Opt] : 
        array_view(object) 
        array operator +( array_view &) 
        array operator -( array_view &) 
        array operator *( array_view &) 
        array operator /( array_view &) 

    cdef cppclass array "triqs::arrays::array" [T,R,Opt] : 
        array()
        array(object)  except +
        array operator +( array_view &) 
        array operator -( array_view &) 
        array operator *( array_view &) 
        array operator /( array_view &) 




