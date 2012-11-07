cdef extern from * :
    cdef cppclass capsule "capsule" [T]:
        capsule(T&)
        object capsule()
        T & content()

cdef class test : 
    cdef test_c _c

    def __init__(self, C_Object) : 
        self.test_c = C_Object


