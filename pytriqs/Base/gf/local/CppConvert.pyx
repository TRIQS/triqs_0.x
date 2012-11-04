#Convert ['triqs::gf::gf_view', 'matrix_view', 'array_view', 'gf_view_freq', 'triqs::arrays::Option::C', 'triqs::arrays::matrix_view']
# TODO: Figure out how many of the pass-by-value copies the compiler can eliminate.

#################### array_view.from_py ####################
{{template_type_declarations}}

cdef extern from "<triqs/arrays.hpp>" namespace "triqs::arrays" : 
    cdef cppclass array_view "array_view" [T,R,Opt] :
        array_view(object)
    ctypedef int TWO   "2" 
    cdef cppclass COrder "C":
        COrder()

@cname("{{cname}}")
cdef array_view[X,TWO,COrder]  {{cname}}(object o) except *:
    print " converting triqs arrays FROM python!"
    return array_view[X,TWO,COrder](o)

#################### C.from_py ####################

cdef extern from "<triqs/arrays.hpp>" namespace "triqs::arrays::Option" : 
    cdef cppclass COrder "triqs::arrays::Option::C":
        COrder()

@cname("{{cname}}")
cdef COrder {{cname}}(object o) except *:
    print " converting triqs arrays FROM python!"
    return COrder()


#################### array_view.to_py ####################

{{template_type_declarations}}

cdef extern from *:
    cdef cppclass array_view "triqs::arrays::array_view" [T,R,Opt] :
        array_view(object)

cdef extern from * namespace triqs::arrays::numpy_interface :
    cdef object array_view_to_python(array_view[T,R,Opt] & )

@cname("{{cname}}")
cdef object {{cname}}(array_view[T,R,Opt] & v):
    print " converting triqs arrays TO python!"
    return array_view_to_python(v)

#################### matrix_view.from_py ####################
{{template_type_declarations}}

cdef extern from "<triqs/arrays.hpp>" namespace "triqs::arrays" : 
    cdef cppclass matrix_view "triqs::arrays::matrix_view" [T,Opt] :
        matrix_view(object) except +
    cdef cppclass COrder "triqs::arrays::Option::C":
        COrder()

@cname("{{cname}}")
cdef matrix_view[X,COrder]  {{cname}}(object o) except *:
    print " converting triqs arrays FROM python!"
    return matrix_view[X,COrder](o)

#################### matrix_view.to_py ####################

{{template_type_declarations}}

cdef extern from *:
    cdef cppclass matrix_view "triqs::arrays::matrix_view" [T,R,Opt] :
        matrix_view(object)

cdef extern from * namespace triqs::arrays::numpy_interface :
    cdef object matrix_view_to_python(matrix_view[T,R,Opt] & )

@cname("{{cname}}")
cdef object {{cname}}(matrix_view[T,R,Opt] & v):
    print " converting triqs arrays TO python!"
    return matrix_view_to_python(v)
#
#
#################### green_fun.from_py ####################

cdef extern from *:
    cdef cppclass green_fun_c "green_fun" :
        green_fun_c(int)
        green_fun_c(green_fun_c)


@cname("{{cname}}")
cdef green_fun_c  {{cname}}(object o) except *:
    print " converting green_fun !"
    return green_fun_c(19) #deref((<green_fun_c?>o).ptr)


#################### gf_view.from_py ####################

cdef extern from * :
    cdef cppclass gf_view_freq_c "triqs::gf::gf_view<triqs::gf::matsubara_freq>" :
        gf_view_freq_c()
        gf_view_freq_c(gf_view_freq_c)

    cdef cppclass gf_view "triqs::gf::gf_view" [D] :
        gf_view()
        gf_view(gf_view)

    cdef cppclass matsubara_freq_desc "triqs::gf::matsubara_freq" : 
        matsubara_freq_desc()

cdef class GFBloc_ImFreq_cython:
    cdef gf_view_freq_c _c
 
#from cython.operator cimport dereference as deref, preincrement as inc

@cname("{{cname}}")
cdef gf_view[matsubara_freq_desc]  {{cname}}(object o) except *:
    print " converting gf_view_freq !"
    #return gf_view_freq(1) #o._get_c()
    return (<GFBloc_ImFreq_cython>o)._c
    return gf_view[D] ( o.mesh, o._data, o._tail) 

#################### gf_view_freq.to_py ####################

cdef extern from * :
    cdef cppclass gf_view_freq "triqs::gf::gf_view<triqs_gf_matsubara_freq_desc>" :
        gf_view_freq()
        gf_view_freq(matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view, nothing) #except +
        void rebind( gf_view_freq&)
        matsubara_freq_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

    cdef gf_view_freq matsubara_freq_make_gf "triqs::gf::matsubara_freq::make_gf" (matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +


@cname("{{cname}}")
cdef object {{cname}}(gf_view_freq & v):
        return []

