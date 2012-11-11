import numpy
from ArrayViewWithIndexConverter import ArrayViewWithIndexConverter,_IndicesConverter
from TailGF import TailGf 
import impl_plot
import impl_lazy

cdef class GfImFreq :
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
        
        c_obj = d.pop('C_Object', None)
        if c_obj :
            assert d == {}
            self._c = extractor [gf_imfreq] (c_obj) () 
            return

        if 'Mesh' not in d : 
            if 'Beta' not in d : raise ValueError, "Beta not provided"
            Beta = float(d.pop('Beta'))
            Nmax = d.pop('NFreqMatsubara',1025)
            stat = d.pop('Statistic','F') # GF_Statistic.Fermion
            sh = 1 if stat== 'F' else 0 # GF_Statistic.Fermion else 0
            d['Mesh'] = MeshImFreq(Beta,'F',Nmax)

        indL = list ( d.pop('IndicesL',()) or d.pop('Indices',()) )
        indR = list ( d.pop('IndicesR',()) or d.pop('Indices',()) )
        assert self.N1 == len(indL) and self.N2 == len(indR)

        cdef MeshImFreq mesh = d.pop('Mesh')
        data_raw = d.pop('Data') if 'Data' in d else numpy.zeros((len(indL),len(indR),len(mesh)), numpy.complex)
        tail= d.pop('Tail') if 'Tail' in d else TailGf(OrderMin=-1, size=10, IndicesL=indL, IndicesR=indR)

        self.Name = d.pop('Name','g')
        
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        self._c =  gf_imfreq ( mesh._c, array_view[dcomplex,THREE,COrder](data_raw), (<TailGf_cython>tail)._c , nothing()) # add here the indices ...
        self._myIndicesGFBlocL = _IndicesConverter(indL)
        self._myIndicesGFBlocR = _IndicesConverter(indR)
        # end of construction ...
    
    # Access to elements of _c, only via C++
    property mesh : 
        """Mesh"""
        def __get__(self): return make_MeshImFreq (self._c.mesh())
    
    property _tail : 
        def __get__(self): return make_TailGf (self._c.singularity_view()) 

    property N1 : 
        def __get__(self): self._c.data_view().dim(0)

    property N2 : 
        def __get__(self): self._c.data_view().dim(1)

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
            return ArrayViewWithIndexConverter(self._c.data_view().to_python(), self._indL, self._indR, None)
        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            a[:,:,:] = value

    property indices : 
        """Indices ..."""
        def __get__(self) : 
            if self._IndicesL != self._IndicesR : raise RuntimeError, "Indices R and L are not the same. I can not give you the Indices"
            for ind in self._IndicesL: 
                yield ind

    property indicesL : 
        """Indices ..."""
        def __get__(self) : 
            for ind in self._IndicesL: 
                yield ind
    
    property indicesR : 
        """Indices ..."""
        def __get__(self) : 
            for ind in self._IndicesR: 
                yield ind

    def _array_with_indices(self) :
        """ 
        Returns ArrayViewWithIndexConverter with : 
            * the indices of the Green function
            * the correct size
            * an array initialized to 0
        """
        return ArrayViewWithIndexConverter(A = numpy.zeros((self.N1,self.N2), dtype = numpy.complex_),
                                           IndicesL = self._IndicesL, IndicesR = self._IndicesR)

    # -------------- Fourier ----------------------------

    def setFromFourierOf(self,GfImTime_cython gt) :
        """Fills self with the Fourier transform of gt"""
        self._c = lazy_fourier( gt._c )

    #-------------   COPY ----------------------------------------

    def copy (self) : 
        r = make_GfImFreq( clone_gf_imfreq(self._c))
        r.Name, r.Note = self.Name, self.Note
        return r
        
    def copy_from(self, GfImFreq G) :
        # Check that dimensions are ok ...
        self._c = G._c

    #---------------------   [  ] operator        ------------------------------------------
    
    def __getitem__(self,key):
        """Key is a tuple of index (n1,n2) as defined at construction"""
        if len(key) !=2 : 
            raise KeyError, "[ ] must be given two arguments"
        try :
            sl1 = self._myIndicesGFBlocL.convertToNumpyIndex(key[0])
        except IndexError, ValueError:
            raise IndexError, "Indices %s incorrect. Indices are %s"%(key[0],list(self._IndicesL))

        try :
            sl2 = self._myIndicesGFBlocR.convertToNumpyIndex(key[1])
        except IndexError, ValueError:
            raise IndexError, "Indices %s incorrect. Indices are %s"%(key[1],list(self._IndicesR))

        return self.__class__(IndicesL = self._IndicesL[sl1],
                              IndicesR = self._IndicesR[sl2],
                              Name = self.Name,
                              Mesh = self._mesh,
                              Data = self._data.array[sl1,sl2,:],
                              Tail = self._tail._make_slice(sl1,sl2))

    def __setitem__(self,key,val):
        g = self.__getitem__(key)
        g <<= val

    #------------- Iteration ------------------------------------

    def __iter__(self) :
        for i in self._IndicesL : 
            for j in self._IndicesR :
                b =self[i,j]
                b.Name = "%s_%s_%s"%(self.Name if hasattr(self,'Name') else '',i,j)
                yield i,j,b
 
    #-------------- Reduction -------------------------------

    #def __reduce__(self):
    #    return call_factory_from_dict, (self.__class__,self.__reduce_to_dict__())
        
    def __reduce_to_dict__(self):
        indL = repr(tuple(self.IndicesL))
        indR = repr(tuple(self.IndicesR))
        assert(eval(indL)==tuple(self.IndicesL))
        assert(eval(indR)==tuple(self.IndicesR))
        return {'IndicesL' : indL,
                'IndicesR' : indR,
                'Data' : self._data,
                'Mesh' : self._mesh,
                'Tail' : self._tail,
                'Name' : self.Name,
                'Note' : self.Note }

    #---------------- Repr, str ---------------------------------

    def __str__ (self) : 
        return self.Name if self.Name else repr(self)
    
    def __repr__(self) : 
        return """%s %s : IndicesL = %s, IndicesR = %s"""%(self.__class__.__name__, self.Name,
          [x for x in self.IndicesL], [x for x in self.IndicesR])

    #--------------   PLOT   ---------------------------------------

    property real : 
        """Use self.real in a plot to plot only the real part"""
        def __get__ (self): return impl_plot.Plot_Wrapper_Partial_Reduce(self,RI='R')

    property imag : 
        """Use self.imag in a plot to plot only the imag part"""
        def __get__ (self): return impl_plot.Plot_Wrapper_Partial_Reduce(self,RI='I')
    
    def _plot_(self, OptionsDict):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        return impl_plot.plot_base( OptionsDict,  r'$\omega_n$', 
                lambda name : r'%s$(i\omega_n)$'%name, True, [x.imag for x in self._mesh] )
    
    #------------------
    
    def x_data_view (self, x_window = None, flatten_y = False) : 
        """
         :param x_window: the window of x variable (omega/omega_n/t/tau) for which data is requested
                          if None, take the full window
         :param flatten_y: If the Green function is of size (1,1) flatten the array as a 1d array
         :rtype: a tuple (X, data) where 
                 * X is a 1d numpy of the x variable inside the window requested
                 * data is a 3d numpy array of dim (:,:, len(X)), the corresponding slice of data
                   If flatten_y is True and dim is (1,1,*), returns a 1d numpy
        """
        return impl_plot.x_data_view(self, [x.imag for x in self._mesh] , x_window, flatten_y)
        
    #--------  LAZY expression system -----------------------------------------
    def __lazy_expr_eval_context__(self) : return impl_lazy.ctx(self)

    #def __le__(self, other) :
    #    """ Forbidden : to avoid typo with <<="""
    #    raise RuntimeError, " Operator <= not defined "
    
    def __ilshift__(self, A): 
        """ A can be two things :
          * G <<= any_GF_Initializers will init the GFBloc with the initializer
          * G <<= g2 where g2 is a GFBloc will copy g2 into self
        """
        return impl_lazy.ilshift_impl(self,A)

    #--------------------  Arithmetic operations  ---------------------------------

    def __iadd__(self, GfImFreq arg):
        self._c += arg._c
        return self

    def __add__(self, GfImFreq y):
        c = self.copy()
        c += y
        return c

    def __isub__(self, GfImFreq arg):
        self._c -= arg._c
        return self

    def __sub__(self,GfImFreq y):
        c = self.copy()
        c -= y
        return c

    def __imul__(self,arg):
        """ If arg is a scalar, simple scalar multiplication
            If arg is a GF (any object with _data and _tail as in GF), they it is a matrix multiplication, slice by slice
        """
        if type(self).__name__ != 'GfImFreq' : self, arg = arg, self
        n = type(arg).__name__
        if n == 'GfImFreq' :
            y = <GfImFreq?>arg
            self._c = self._c * as_gf_imfreq(y)
        elif n in ['float','int', 'complex'] : 
            self._c = <dcomplex>arg * self._c
        else : 
            a = numpy.array(arg, numpy.complex)
            self._c = self._c * matrix[dcomplex,COrder](a) 
        return self

    def __mul__(self,y):
        
        # This is for time only ?? What about real freq, im freq ??
        #if hasattr(y,"_data") :
        #    c = self.copy_with_new_stat(GF_Statistic.Boson if self._mesh.Statistic == y._mesh.Statistic else GF_Statistic.Fermion)
        #else:
        c = self.copy()
        c *= y
        return c

    def __idiv__(self,arg):
        if type(self).__name__ != 'GfImFreq' : 
            raise RuntimeError, "Can not divide by an GfImFreq"
        if type(arg).__name__  in ['float','int', 'complex'] : 
            self._c = self._c / <dcomplex>arg
        else : 
            raise RuntimeError, " argument type not recognized for %s"%arg
        return self

    def __div__(self,y):
        c = self.copy()
        c /= y
        return c

    #--------------   OTHER OPERATIONS -----------------------------------------------------

    def from_L_G_R (self, L,G,R):
        """ For all argument, replace the matrix by L *matrix * R"""
        self._c = matrix[dcomplex,COrder](L) * self._c * matrix[dcomplex,COrder](R) 

    def invert(self) : 
        """Invert the matrix for all arguments"""
        self._c = inverse (self._c)

    def replace_by_tail(self,start) : 
        d = self._data_raw
        t = self._tail
        for n, om in enumerate(self.mesh) : 
            if n >= start : d[:,:,n] = t(om).array

    def transpose(self):
        """Transposes the GF Bloc: return a new transposed view"""
        r = make_GfImFreq ( transpose (self._c) )
        self.Name=self.Name + '(t)'
        return r

    def conjugate(self):
        """Complex conjugate of the GF Bloc. It follow the policy of numpy and
        make a copy only if the Green function is complex valued"""
        r = make_GfImFreq ( conjugate (self._c) )
        self.Name=self.Name + '*'
        return r
    
    #------------------  Density -----------------------------------
    
    def total_density(self) :
        """Trace density"""
        return numpy.trace(self.density())

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (GfImFreq)

