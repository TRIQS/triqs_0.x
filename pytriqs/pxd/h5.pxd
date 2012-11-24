from libcpp.string cimport string as std_string
#cdef extern from "<triqs/utility/h5_tools.hpp>" namespace "triqs::h5" : 
cdef extern from "<triqs/arrays/h5/group_or_file.hpp>" : #namespace "triqs::arrays::h5" : 

    cdef cppclass h5_group_or_file "triqs::arrays::h5::group_or_file" : 
        h5_group_or_file (int, bint)
        h5_group_or_file ()

    cdef cppclass h5_extractor "triqs::arrays::h5::h5_extractor" [T] : 
        h5_extractor()
        T & operator()( h5_group_or_file &, std_string)

cdef inline h5_group_or_file make_h5_group_or_file (gr) :
        import h5py
        return h5_group_or_file(gr.id.id, type(gr) == h5py._hl.group.Group)
