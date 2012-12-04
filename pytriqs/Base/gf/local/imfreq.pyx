# ----------- Mesh  --------------------------
cdef class MeshImFreq: 
    cdef mesh_imfreq  _c

    def __init__(self, Beta, stat, int Nmax): 
        self._c =  make_mesh_imfreq(Beta,{'F' :Fermion, 'B' : Boson}[stat] ,Nmax) 
    
    def __len__ (self) : return self._c.size()
    
    property beta : 
        """Inverse temperature"""
        def __get__(self): return self._c.domain().beta
    
    property statistic : 
        def __get__(self): return 'F' if self._c.domain().statistic==Fermion else 'B'
    
    def __iter__(self) : # I use the C++ generator !
        cdef mesh_pt_generator[mesh_imfreq ] g = mesh_pt_generator[mesh_imfreq ](&self._c)
        while not g.at_end() : 
            yield g.to_point()
            g.increment()

    def __richcmp__(MeshImFreq self, MeshImFreq other,int op) : 
        if op ==2 : # ==
            return self._c == other._c

# C -> Python 
cdef inline make_MeshImFreq ( mesh_imfreq x) :
    return MeshImFreq( x.domain().beta, 'F', x.size() )
    #return MeshImFreq(encapsulated_c_object = encapsulate (&x))

# ----------- GF --------------------------

cdef class GfImFreq(_ImplGfLocal) :
    cdef gf_imfreq _c
    def __init__(self, **d):
        """
        The constructor have two variants : you can either provide the mesh in
        Matsubara frequencies yourself, or give the parameters to build it.
        All parameters must be given with keyword arguments.

        GfImFreq(Indices, Beta, Statistic, NFreqMatsubara,  Data, Tail, Name)

               * ``Indices``:  a list of indices names of the block
               * ``Beta``:  Inverse Temperature 
               * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
               * ``NFreqMatsubara``:  Number of Matsubara frequencies
               * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NFreqMatsubara) representing the value of the Green function on the mesh. 
               * ``Tail``:  the tail 
               * ``Name``:  a name of the GF

        If you already have the mesh, you can use a simpler version :

        GfImFreq(Indices, Mesh, Data, Tail, Name)
            
               * ``Indices``:  a list of indices names of the block
               * ``Mesh``:  a MeshGF object, such that Mesh.TypeGF== GF_Type.Imaginary_Frequency 
               * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NFreqMatsubara) representing the value of the Green function on the mesh. 
               * ``Tail``:  the tail 
               * ``Name``:  a name of the GF

        .. warning::
        The Green function take a **view** of the array Data, and a **reference** to the Tail.
        """
        c_obj = d.pop('encapsulated_c_object', None)
        if c_obj :
            assert d == {}, "Internal error : encapsulated_c_object must be the only argument"
            self._c = extractor [gf_imfreq] (c_obj) ()
            self._IndicesL = self._c.indices()[0]
            self._IndicesR = self._c.indices()[1]
            return 

        bss = d.pop('boost_serialization_string', None)
        if bss :
            assert d == {}, "Internal error : boost_serialization_string must be the only argument"
            boost_unserialize_into(<std_string>bss,self._c) 
            return 

        _ImplGfLocal.__init__(self, d)

        cdef MeshImFreq mesh = d.pop('Mesh',None)
        if mesh is None : # 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NFreqMatsubara',1025)
            stat = d.pop('Statistic','F') 
            sh = 1 if stat== 'F' else 0 
            #d['Mesh'] = MeshImFreq(Beta,'F',Nmax)
            mesh = MeshImFreq(Beta,'F',Nmax)

        self.dtype = numpy.complex_
        data_raw = d.pop('Data') if 'Data' in d else numpy.zeros((len(self._IndicesL),len(self._IndicesR),len(mesh)), self.dtype )
        
        cdef TailGf tail= d.pop('Tail') if 'Tail' in d else TailGf(OrderMin=-1, size=10, IndicesL=self._IndicesL, IndicesR=self._IndicesR)

        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        
        self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data_raw), tail._c , nothing(), make_c_indices(self) ) 
       
        # object to keep that will never change
        #self.mesh = make_MeshImFreq (self._c.mesh())
        #self.N1, self.N2 = self._c.data_view().shape(0), self._c.data_view().shape(1)
        # end of construction ...
    
    # Access to elements of _c, only via C++
    property mesh : 
        """Mesh"""
        def __get__(self): return make_MeshImFreq (self._c.mesh())
    
    property tail : 
        def __get__(self): return make_TailGf (self._c.singularity_view()) 
        def __set__(self,TailGf t): 
            assert (self.N1, self.N2, self._c.singularity_view().size()) == (t.N1, t.N2, t.size)
            cdef tail t2 = self._c.singularity_view()
            t2 = t._c 

    property N1 : 
        def __get__(self): return self._c.data_view().shape(0)

    property N2 : 
        def __get__(self): return self._c.data_view().shape(1)

    property data : 
        """Access to the data array"""
        def __get__(self) : 
            return self._c.data_view().to_python()
        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            a[:,:,:] = value
   
    def __typename(self) : return "GfImFreq"

    #-------------- Reduction -------------------------------

    def __reduce__(self):
        return py_deserialize, (self.__class__,boost_serialize(self._c),)

    # -------------- HDF5 ----------------------------

    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group_or_file(gr), key, self._c)

    # -------------- Fourier ----------------------------

    def set_from_fourier_of(self,GfImTime gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( gt._c )

    #-------------   COPY ----------------------------------------

    def copy (self) : 
        r = make_GfImFreq( clone_gf_imfreq(self._c))
        r.Name = self.Name
        return r
        
    def copy_from(self, GfImFreq G) :
        # Check that dimensions are ok ...
        self._c = G._c

    #--------------   PLOT   ---------------------------------------
   
    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        return self._plot_base( OptionsDict,  r'$\omega_n$', 
                lambda name : r'%s$(i\omega_n)$'%name, True, [x.imag for x in self.mesh] )
    #--------------   OTHER OPERATIONS -----------------------------------------------------

    def invert(self) : 
        """Invert the matrix for all arguments"""
        pass
        #self._c = inverse_c (self._c)

    def replace_by_tail(self,start) : 
        d = self.data
        t = self.tail
        for n, om in enumerate(self.mesh) : 
            if n >= start : d[:,:,n] = t(om).array

#----------------  Reading from h5 ---------------------------------------

def h5_read_GfImFreq ( gr, std_string key) : 
    return make_GfImFreq( h5_extractor[gf_imfreq]()(make_h5_group_or_file(gr),key))

#----------------  Convertions functions ---------------------------------------

# Python -> C
cdef gf_imfreq  as_gf_imfreq (g) except +: 
    return (<GfImFreq?>g)._c

# C -> Python. Do NOT add except +
cdef make_GfImFreq ( gf_imfreq  x)  :
        return GfImFreq(encapsulated_c_object = encapsulate (&x))

# Python -> C for blocks
cdef gf_block_imfreq  as_gf_block_imfreq (G) except +:
        cdef vector[gf_imfreq] v_c
        for item in G:
            v_c.push_back(as_gf_imfreq(item))
        return make_gf_block_imfreq (v_c)

# C -> Python for block
cdef make_BlockGfImFreq (gf_block_imfreq G) :
    gl = []
    name_list = G.mesh().domain().names()
    cdef int i =0
    for n in name_list:
        gl.append( make_GfImFreq(G[i] ) )
    return GF( NameList = name_list, BlockList = gl)

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (GfImFreq, read_fun = h5_read_GfImFreq)

