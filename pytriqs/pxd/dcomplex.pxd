cdef extern from "<complex>" namespace "std":
    cdef cppclass dcomplex "std::complex<double>":
         dcomplex()
         dcomplex(dcomplex &)
         dcomplex(double,double)
         double real, imag

# Python -> C
cdef inline dcomplex as_dcomplex (a) : 
    x = complex(a)
    return dcomplex(a.real, a.imag)

# C -> Python 
cdef inline make_dcomplex (dcomplex z) :
    return complex(z.real, z.imag)





