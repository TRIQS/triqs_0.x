# ----------- Domain  --------------------------
#cdef class DomainMatsubaraTime:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshMatsubaraTime: 
    #cdef matsubara_time_mesh_c _c
    
    def __init__(self, Beta, stat, int Nmax): 
        self._c = matsubara_time_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) :
        cdef mesh_pt_generator[matsubara_time_mesh_c] g = mesh_pt_generator[matsubara_time_mesh_c](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- The GF  --------------------------

cdef class GFBloc_ImTime_cython:
    #cdef gf_im_time_c _c
    #cdef object _mesh
    def __init__(self, MeshMatsubaraTime mesh, data, TailGF_c tail):
        self._c =  gf_im_time_c ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c, nothing() )
        self._mesh = mesh
    
    def setFromInverseFourierOf(self,  gw) : 
        """Fills self with the Inverse Fourier transform of gw"""        
        self._c = lazy_inverse_fourier( (<GFBloc_ImFreq_cython>gw)._c)

