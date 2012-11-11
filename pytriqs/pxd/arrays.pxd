cdef extern from "<triqs/arrays.hpp>" namespace "triqs::arrays::Option" : 
    #ctypedef C COrder 
    #ctypedef Fortran FortranOrder 
    cdef cppclass COrder "triqs::arrays::Option::C":
        COrder()
    cdef cppclass FortranOrder "triqs::arrays::Option::Fortran" :
        FortranOrder()

    
cdef extern from "<triqs/arrays.hpp>" namespace "triqs::arrays" : 
    ctypedef int ONE   "1" 
    ctypedef int TWO   "2" 
    ctypedef int THREE "3" 
    ctypedef int FOUR  "4" 
    ctypedef int FIVE  "5" 
    ctypedef int SIX   "6" 

    #cdef cppclass COrder "triqs::arrays::Option::C":
    #    COrder()
    #cdef cppclass FortranOrder "triqs::arrays::Option::Fortran" :
    #    FortranOrder()

    cdef cppclass array_view "triqs::arrays::array_view" [T,R,Opt] : 
        array_view() 
        array_view(object) except +
        array_view operator +( array_view &) 
        array_view operator -( array_view &) 
        array_view operator *( array_view &) 
        array_view operator /( array_view &) 
        object to_python()
        int shape(int)

    cdef cppclass array "triqs::arrays::array" [T,R,Opt] : 
        array()
        array(object)  except +
        array operator +( array_view &) 
        array operator -( array_view &) 
        array operator *( array_view &) 
        array operator /( array_view &) 
        object to_python()

    cdef cppclass matrix "triqs::arrays::matrix" [T,Opt] :
        matrix()
        matrix(object)  except +
        matrix operator +( matrix_view &)
        matrix operator -( matrix_view &)
        matrix operator /( matrix_view &)
        object to_python()
        
    cdef matrix operator *( matrix &, matrix_view &)
    
    cdef cppclass matrix "triqs::arrays::matrix" [T,Opt] :
        matrix() 
        matrix(object)  except +
 
    cdef cppclass matrix_view "triqs::arrays::matrix_view" [T,Opt] :
        matrix_view() 
        matrix_view(matrix) 
        matrix_view(object)  except +
        #matrix_view( matrix[T,Opt] &)
        object to_python()

    cdef matrix_view operator +( matrix_view & , matrix_view &)
    cdef matrix_view operator -( matrix_view & , matrix_view &)
    cdef matrix_view operator *( matrix_view & , matrix_view &)
    cdef matrix_view operator /( matrix_view & , matrix_view &)
 
    cdef cppclass tqa_vector_view "triqs::arrays::vector_view" [T,R,Opt] :
        tqa_vector_view() 
        tqa_vector_view(object) except +
        tqa_vector_view operator +(tqa_vector_view &)
        tqa_vector_view operator -(tqa_vector_view &)
        tqa_vector_view operator *(tqa_vector_view &)
        tqa_vector_view operator /(tqa_vector_view &)
        object to_python()

    cdef cppclass tqa_vector "triqs::arrays::vector" [T,R,Opt] :
        tqa_vector()
        tqa_vector(object)  except +
        tqa_vector operator +( tqa_vector_view &)
        tqa_vector operator -( tqa_vector_view &)
        tqa_vector operator *( tqa_vector_view &)
        tqa_vector operator /( tqa_vector_view &)
        object to_python()

