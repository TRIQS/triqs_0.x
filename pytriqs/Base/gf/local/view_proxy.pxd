cdef extern from "<triqs/utility/view_tools.hpp>" namespace "triqs" : 
   cdef cppclass view_proxy [GF]:
        view_proxy()
        view_proxy(GF) except +
        void rebind(GF &)
        void operator << (GF &)
        GF & operator ()()


