cdef extern from "<triqs/arrays/h5/group_or_file.hpp>" : #namespace "triqs::arrays::h5" : 

    cdef cppclass h5_group_or_file "triqs::arrays::h5::group_or_file" : 
        #h5_group_or_file (char * filename, int flag)
        h5_group_or_file (int)

#cdef extern int H5F_ACC_TRUNC 
