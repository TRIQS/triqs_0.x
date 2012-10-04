cdef extern from * : 
    cdef cppclass shared_ptr "std::shared_ptr" [T] : 
        shared_ptr(T*)
        shared_ptr()
        T* get()

    ctypedef int ONE   "1" 
    ctypedef int TWO   "2" 
    ctypedef int THREE "3" 
    ctypedef int FOUR  "4" 

    cdef cppclass array_view "triqs::arrays::array_view" [T,R] : 
        array_view(object) #except +
        array operator +( array_view &) 
    
    cdef cppclass array "triqs::arrays::array" [T,R] : 
        array()
        array(object)  except +
        array operator +( array_view &) 




