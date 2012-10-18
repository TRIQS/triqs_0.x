# ##############################################################
#               Matsubara Freq : local GF 
# ##############################################################

# ----------- Descriptor  --------------------------
cdef extern from "triqs/gf/descriptors/matsubara_freq.hpp" namespace "triqs::gf" : 
  
    cdef cppclass matsubara_freq_domain :
        double beta
        statistic_enum statistic
        matsubara_freq_domain ()

    cdef cppclass matsubara_freq_mesh " triqs::gf::matsubara_freq::mesh_t"  :
        double beta
        statistic_enum statistic
        matsubara_freq_mesh ()
        matsubara_freq_mesh (double Beta, statistic_enum s, int Nmax) 
        matsubara_freq_domain & domain()
        long size()
 
# ----------- Domain  --------------------------
cdef class DomainMatsubaraFrequency:
    pass

# ----------- Mesh  --------------------------
cdef class MeshMatsubaraFrequency: 
    cdef matsubara_freq_mesh _c
    def __init__(self, Beta, stat, int Nmax): 
        self._c = matsubara_freq_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    def __len__ (self) : return self._c.size()
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[matsubara_freq_mesh] g = mesh_pt_generator[matsubara_freq_mesh](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()


# ----------- The GF  --------------------------

cdef extern from "triqs/gf/gf.hpp" namespace "triqs::gf" : 
    cdef cppclass gf_view_freq "triqs::gf::gf_view<triqs::gf::matsubara_freq>"  : 
        gf_view_freq()
        gf_view_freq(matsubara_freq_mesh, array_view[dcomplex, THREE,COrder], tail_view) except +
        matsubara_freq_mesh mesh() 
        array_view[dcomplex, THREE,COrder] data_view()
        tail_view auxiliary_data_view() 

cdef extern from "triqs/gf/local/fourier_matsubara.hpp" : 
    gf_view_freq fourier_direct (gf_view_time & )

cdef class GFBloc_ImFreq_c:
    cdef gf_view_freq _c
    cdef object _mesh
    def __init__(self, **args) :
        #if not args : self._c = gf_view_freq () # no argument, empty view
        #else : 
            def deleg (MeshMatsubaraFrequency mesh, data, TailGF_c tail):
                #self._c = gf_view_freq( mesh, data, tail)
                self._c = gf_view_freq ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c )
            deleg(**args) 
    def setFromFourierOf(self, GFBloc_ImTime_c gt) : 
        self._c = fourier_direct( gt._c)

#def fourier_direct3(GFBloc_ImTime gt) : 
#    return fourier_direct(gt)

#def fourier_direct2(GFBloc_ImTime_c gt) : 
#    res = GFBloc_ImFreq_c()
#    res._c = fourier_direct(gt._c)
#   return res

"""
cdef class GFBloc_ImFreq_c:
    cdef gf_view_freq _c
    def __init__(self, **args) :
        if args : self.__init_deleg__(**args)
        else : self._c = gf_view_freq ()
    def __init_deleg__(self, MeshMatsubaraFrequency mesh, data, TailGF_c tail):
        self._c = gf_view_freq ( deref(mesh._c), array_view[dcomplex,THREE,COrder](data), deref(tail._c) )
        #self._c = new gf_view_freq ( deref(mesh._c), array_view[dcomplex,THREE,COrder](data), deref(tail._c) )
    #def __dealloc__(self): del self._c
    def setFromFourierOf(self, GFBloc_ImTime_c gt) : 
        #cdef view_proxy[gf_view_freq] P
        #P.assign( fourier_direct (deref(gt._c)))
        #view_proxy[gf_view_freq]().assign( fourier_direct (deref(gt._c)))
        #view_proxy[gf_view_freq](self._c, fourier_direct (deref(gt._c)))
        #self._c.assign_from (fourier_direct (deref(gt._c)))
        self._c = fourier_direct( gt._c)
        #deref(self._c) = fourier_direct (deref(gt._c))
        #deref(self._c) = fourier_direct (deref(gt._c))

def fourier_direct2(GFBloc_ImTime_c gt) : 
    #cdef gf_view_freq res = fourier_direct(gt._c)
    res = GFBloc_ImFreq_c()
    res._c = fourier_direct(gt._c)
    return res
    #return GFBloc_ImFreq_c( res.mesh(),res.data_view(), res.auxiliary_data_view()) 
"""



