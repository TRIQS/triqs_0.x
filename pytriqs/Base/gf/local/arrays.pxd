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

    cdef cppclass matrix_view "triqs::arrays::matrix_view" [T,R,Opt] :
        matrix_view(object)
        matrix operator +( matrix_view &)
        matrix operator -( matrix_view &)
        matrix operator *( matrix_view &)
        matrix operator /( matrix_view &)

    cdef cppclass matrix "triqs::arrays::matrix" [T,R,Opt] :
        matrix()
        matrix(object)  except +
        matrix operator +( matrix_view &)
        matrix operator -( matrix_view &)
        matrix operator *( matrix_view &)
        matrix operator /( matrix_view &)

    cdef cppclass vector_view "triqs::arrays::vector_view" [T,R,Opt] :
        vector_view(object)
        vector operator +(vector_view &)
        vector operator -(vector_view &)
        vector operator *(vector_view &)
        vector operator /(vector_view &)

    cdef cppclass vector "triqs::arrays::vector" [T,R,Opt] :
        vector()
        vector(object)  except +
        vector operator +( vector_view &)
        vector operator -( vector_view &)
        vector operator *( vector_view &)
        vector operator /( vector_view &)

