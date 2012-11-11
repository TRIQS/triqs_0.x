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
    #return MeshImFreq( x.domain().beta, x.domain().statistic, x.size() )
    #return MeshImFreq(C_Object = encapsulate (&x))

# ----------- GF --------------------------

cdef class GfImFreq(_ImplGfLocal) :
    #cdef gf_imfreq _c
    #object _myIndicesGFBlocL, _myIndicesGFBlocR, Name, dtype

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
        
        c_obj = d.pop('C_Object', None)
        if c_obj :
            assert d == {}
            self._c = extractor [gf_imfreq] (c_obj) () 
            n0,n1 = self._c.data_view().shape(0),self._c.data_view().shape(1)
            _ImplGfLocal.__init__(self, IndicesL = range(n0), IndicesR = range(n1), Name = "")
            return

        if 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NFreqMatsubara',1025)
            stat = d.pop('Statistic','F') 
            sh = 1 if stat== 'F' else 0 
            d['Mesh'] = MeshImFreq(Beta,'F',Nmax)

        self.dtype = numpy.complex_
        indL = list ( d.pop('IndicesL',()) or d.pop('Indices',()) )
        indR = list ( d.pop('IndicesR',()) or indL )

        cdef MeshImFreq mesh = d.pop('Mesh')
        data_raw = d.pop('Data') if 'Data' in d else numpy.zeros((len(indL),len(indR),len(mesh)), self.dtype )
        cdef TailGf tail= d.pop('Tail') if 'Tail' in d else TailGf(OrderMin=-1, size=10, IndicesL=indL, IndicesR=indR)

        _ImplGfLocal.__init__(self, IndicesL = indL, IndicesR = indR, Name =  d.pop('Name','g'))
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        
        self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data_raw), tail._c , nothing()) # add here the indices ...
        assert self.N1 == len(indL) and self.N2 == len(indR)
        # end of construction ...
    
    # Access to elements of _c, only via C++
    property mesh : 
        """Mesh"""
        def __get__(self): return make_MeshImFreq (self._c.mesh())
    
    property _tail : 
        def __get__(self): return make_TailGf (self._c.singularity_view()) 
        def __set__(self,TailGf t): 
                assert (self.N1, self.N2, self._c.singularity_view().size()) == (t.N1, t.N2, t.size)
                cdef tail t2 = self._c.singularity_view()
                t2 = t._c 

    property N1 : 
        def __get__(self): return self._c.data_view().shape(0)

    property N2 : 
        def __get__(self): return self._c.data_view().shape(1)

    property _data_raw : 
        """Access to the data array"""
        def __get__(self) : 
            return self._c.data_view().to_python()
        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            a[:,:,:] = value
    
    property _data : 
        """Access to the data array"""
        def __get__(self) :
            return ArrayViewWithIndexConverter(self._c.data_view().to_python(), self.indicesL, self.indicesR, None)
        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            a[:,:,:] = value

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
                lambda name : r'%s$(i\omega_n)$'%name, True, [x.imag for x in self._mesh] )
    
    #--------------------  Arithmetic operations  ---------------------------------

    def __iadd__(self, GfImFreq arg):
        self._c = self._c + arg._c
        return self

    def __add__(self, GfImFreq y):
        c = self.copy()
        c += y
        return c

    def __isub__(self, GfImFreq arg):
        self._c = self._c - arg._c
        return self

    def __sub__(self,GfImFreq y):
        c = self.copy()
        c -= y
        return c

    def __imul__(self,arg):
        """ If arg is a scalar, simple scalar multiplication
            If arg is a GF (any object with _data and _tail as in GF), they it is a matrix multiplication, slice by slice
        """
        n = type(arg).__name__
        if n == 'GfImFreq' :
            self._c = self._c * (<GfImFreq?>arg)._c
        elif n in ['float','int', 'complex'] : 
            self._c = as_dcomplex(arg) * self._c
        else : 
            raise RuntimeError, " argument type not recognized in imul for %s"%arg
        return self

    def __mul_impl__(self, arg, s) : 
        cdef GfImFreq res = self.copy()
        n = type(arg).__name__
        cdef matrix_view [dcomplex,COrder] a 
        if n == 'GfImFreq' :
            res._c =  self._c * (<GfImFreq?>arg)._c
        elif n in ['float','int', 'complex'] : 
            res._c = as_dcomplex(arg) * self._c
        else : 
            a= matrix_view[dcomplex,COrder](matrix[dcomplex,COrder](numpy.array(arg, self.dtype)))
            #res._c =  a * self._c  if s else self._c *a
        return res

    def __mul__(self,arg):
        s = type(self).__name__ != 'GfImFreq' 
        return self.__mul_impl__(self, arg, s) if not s else self.__mul_impl__(self, arg, s)

    def __idiv__(self,arg):
        cdef GfImFreq me = self
        me._c = me._c / as_dcomplex(arg)
        return self

    def __div_impl_(self, arg, s):
        if s : raise RuntimeError, "Can not divide by an GfImFreq"
        cdef GfImFreq res = self.copy()
        if type(arg).__name__  in ['float','int', 'complex'] : 
            res._c = self._c / as_dcomplex(arg)
        else : 
            raise RuntimeError, " argument type not recognized for %s"%arg
        return res

    def __div__(self,arg):
        s = type(self).__name__ != 'GfImFreq' 
        return self.__div_impl__(self, arg, s) if not s else self.__div_impl__(self, arg, s)

    #--------------   OTHER OPERATIONS -----------------------------------------------------

    def from_L_G_R (self, L,G,R):
        """ For all argument, replace the matrix by L *matrix * R"""
        pass
        #self._c = matrix_view[dcomplex,COrder](L) * self._c * matrix_view[dcomplex,COrder](R) 

    def invert(self) : 
        """Invert the matrix for all arguments"""
        pass
        #self._c = inverse (self._c)

    def replace_by_tail(self,start) : 
        d = self._data_raw
        t = self._tail
        for n, om in enumerate(self.mesh) : 
            if n >= start : d[:,:,n] = t(om).array

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (GfImFreq)

