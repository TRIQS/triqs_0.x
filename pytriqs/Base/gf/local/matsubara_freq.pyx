# ----------- Domain  --------------------------
#cdef class DomainMatsubaraFrequency:
#    pass

# ----------- Mesh  --------------------------
cdef class MeshMatsubaraFrequency: 
    #cdef matsubara_freq_mesh_c  _c

    def __init__(self, Beta, stat, int Nmax): 
        self._c =  matsubara_freq_make_mesh(Beta,{ 'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property Beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property Statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[matsubara_freq_mesh_c ] g = mesh_pt_generator[matsubara_freq_mesh_c ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

# ----------- GF  --------------------------


cdef extern from * :
    object PyCObject_FromVoidPtr(void* cobj, void *)
    void* PyCObject_AsVoidPtr(object self)

from GFBloc_ImFreq import GFBloc_ImFreq 

cdef inline toPy_gf_im_freq_c ( gf_im_freq_c & c) : 
    return GFBloc_ImFreq( C_Object = PyCObject_FromVoidPtr (&c, NULL))
    #return GFBloc_ImFreq_cython( C_Object = c.address_as_opaque_python_object())

cdef class GFBloc_ImFreq_cython:
    #cdef object _mesh
    #cdef gf_im_freq_c _c
    
    def __init__(self, *l, **d):
        def  deleg (MeshMatsubaraFrequency mesh, data, TailGF_c tail):
            self._c =  gf_im_freq_c ( mesh._c, array_view[dcomplex,THREE,COrder](data), tail._c , nothing())
            self._mesh = mesh
        if len(l)==0 and d.keys()==['C_Object'] : 
            self._c = deref(<gf_im_freq_c *>PyCObject_AsVoidPtr(d['C_Object']))
            self._mesh = mesh
            #self._c = gf_im_freq_c ( d['C_Object'],0)
        else : 
            deleg(*l,**d)

    def test(self) :
        return toPy_gf_im_freq_c(self._c)
 
    #def pr(self) :
    #    print &self._c

    def setFromFourierOf(self, gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( (<GFBloc_ImTime_cython >gt)._c )

# ----------- Convertion function --------------------------


