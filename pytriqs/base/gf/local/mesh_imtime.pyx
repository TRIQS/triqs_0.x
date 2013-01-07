cdef class MeshImTime: 
    cdef mesh_imtime  _c

    def __init__(self, Beta, stat, int n_max): 
        self._c =  make_mesh_imtime(Beta,{'F' :Fermion, 'B' : Boson}[stat] ,n_max) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[mesh_imtime ] g = mesh_pt_generator[mesh_imtime ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

    def __richcmp__(MeshImTime self, MeshImTime other,int op) : 
        if op ==2 : # ==
            return self._c == other._c

    def __reduce__(self):
        return self.__class__, (self.beta, self.statistic, len(self))

# C -> Python 
cdef inline make_MeshImTime ( mesh_imtime x) :
    return MeshImTime( x.domain().beta, 'F', x.size() )

