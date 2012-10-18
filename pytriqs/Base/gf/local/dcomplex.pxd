cdef extern from "<complex>" namespace "std":
    cdef cppclass dcomplex "std::complex<double>":
         complex(double,double)


