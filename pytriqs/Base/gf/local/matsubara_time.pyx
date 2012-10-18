# ##############################################################
#               Matsubara Time : local GF 
# ##############################################################

# ----------- Descriptor  --------------------------
cdef extern from "triqs/gf/descriptors/matsubara_time.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_time_domain :
        double beta
        statistic_enum statistic
        matsubara_time_domain ()

    cdef cppclass matsubara_time_mesh " triqs::gf::matsubara_time::mesh_t"  :
        double beta
        statistic_enum statistic
        matsubara_time_mesh () 
        matsubara_time_mesh (double Beta, statistic_enum s, int Nmax) 
        matsubara_time_domain & domain()
        long size()

# ----------- Domain  --------------------------
cdef class DomainMatsubaraTime:
    pass

## ----------- Mesh  --------------------------
cdef class MeshMatsubaraTime: 
    cdef matsubara_time_mesh _c
    def __init__(self, Beta, stat, int Nmax): 
        self._c = matsubara_time_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    def __len__ (self) : return self._c.size()
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    def __iter__(self) :
        cdef mesh_pt_generator[matsubara_time_mesh] g = mesh_pt_generator[matsubara_time_mesh](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- The GF  --------------------------

cdef extern from "triqs/gf/gf.hpp" namespace "triqs::gf" : 
    cdef cppclass gf_view_time "triqs::gf::gf_view<triqs::gf::matsubara_time>"  : 
        gf_view_time()
        gf_view_time(matsubara_time_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +

cdef class GFBloc_ImTime_c:
    cdef gf_view_time _c
    def __init__(self, MeshMatsubaraTime mesh, data, TailGF_c tail):
        self._c = gf_view_time ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c )
