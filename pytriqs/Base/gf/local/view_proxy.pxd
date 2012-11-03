cdef extern from "<triqs/utility/view_tools.hpp>" namespace "triqs" : 
   cdef cppclass view_proxy [GF]:
        view_proxy()
        void rebind(GF &)
        void operator << (GF &)
        GF & operator ()()


