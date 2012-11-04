# ----------- Domain  --------------------------
#cdef class DomainMatsubaraFrequency:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshMatsubaraFrequency: 
    cdef matsubara_freq_mesh  _c

    def __init__(self, Beta, stat, int Nmax): 
        self._c =  matsubara_freq_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[matsubara_freq_mesh ] g = mesh_pt_generator[matsubara_freq_mesh ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- GF  --------------------------

cdef class GFBloc_ImFreq_cython:
    cdef object _mesh
    cdef gf_view_freq _c

    def __init__(self, MeshMatsubaraFrequency mesh, data, TailGF_c tail):
            self._c =  gf_view_freq ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing())
            self._mesh = mesh

    def setFromFourierOf(self, GFBloc_ImTime_cython gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( gt._c )


