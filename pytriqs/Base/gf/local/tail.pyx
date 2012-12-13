cdef class TailGf:
    cdef tail _c
    def __init__(self, **d):
        """
        TailGf ( data, order_min)
        TailGf ( N1, N2, size, order_min) 

              * ``Indices``:  a list of indices names of the block
              * ``Beta``:  Inverse Temperature 
              * ``Statistic``:  GF_Statistic.Fermion [default] or GF_Statistic.Boson
              * ``NFreqMatsubara``:  Number of Matsubara frequencies
              * ``Data``:   A numpy array of dimensions (len(Indices),len(Indices),NFreqMatsubara) representing the value of the Green function on the mesh. 
              * ``Tail``:  the tail 
              * ``Name``:  a name of the GF

        """
        a = d.pop('data',None)
        if a==None : 
            (N1, N2), s = d.pop('shape'), d.pop('size')
            a = numpy.zeros((N1,N2,s) ,numpy.complex) 
        omin = d.pop('order_min')
        assert len(d) ==0, "Unknown parameters in TailGf constructions %s"%d.keys() 
        self._c =  tail( array_view[dcomplex,THREE,COrder](a), omin)
    
    property data : 
        """Access to the data array"""
        def __get__(self) : 
            return self._c.data_view().to_python()

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            if hasattr(value,'shape') : 
                if a.shape[:2] != value.shape[:2] : 
                    raise RuntimeError, "shape mismatch"
                m = min(value.shape[2],a.shape[2]) 
                if a.shape[2] > value.shape[2] : 
                    a[:,:,m:] =0
                a[:,:, 0:m] = value[:,:, 0:m]
            else : 
                a[...] = value

    property shape : 
        def __get__(self) : return self.data.shape[:2] 

    property order_min : 
        """Min order of the expansion"""
        def __get__(self) : return self._c.order_min()

    property order_max : 
        """Max order of the expansion"""
        def __get__(self) : return self._c.order_max()

    property N1 : 
        def __get__(self): return self._c.data_view().shape(0)

    property N2 : 
        def __get__(self): return self._c.data_view().shape(1)

    property size : 
        """Length of the expansion"""
        def __get__(self) : return self._c.size()
    
    def copy(self) : 
        return self.__class__(data = self.data.copy(), order_min = self.order_min)

    # Should I do more compatibility checks?
    def copy_from(self, TailGf T) :
        self._c = T._c

    def _make_slice(self, sl1, sl2):
        return self.__class__(data = self.data[sl1,sl2,:], order_min = self.order_min)

    def __repr__ (self) :
        return string.join([ "%s"%self[r]+ (" /" if r<0 else "") + " Om^%s"%(abs(r)) for r in range(self.order_min, self.order_max+1) ] , " + ")

    def __getitem__(self,i) :
        """Returns the i-th coefficient of the expansion, or order Om^i"""
        if not self.has_coef(i) : raise IndexError, "Index %s is out of range"%i
        return self.data[:,:,i-self.order_min]

    def __setitem__(self,i, val) :
        """Sets the i-th coefficient of the expansion, or order Om^i"""
        if not self.has_coef(i) : raise IndexError, "Index %s is out of range"%i
        self.data[:,:,i-self.order_min] = val

    def has_coef(self, i):
        return (i >= self.order_min) and (i <= self.order_max)

    def __call__(self, n) :
        if not self.has_coef(n) : raise IndexError, "Index %s is out of range"%n
        return self._c(n).to_python()   

    def __reduce__(self):
        return (lambda cls, d : cls(**d)) , (self.__class__,self.__reduce_to_dict__())
  
    def invert(self) :
        self._c = 1.0/self._c
    
    #########      arithmetic operations    #################

    def __iadd__(self, TailGf arg):
        self._c = self._c + arg._c
        return self

    def __add__(self, TailGf y):
        c = self.copy()
        c += y
        return c

    def __isub__(self, TailGf arg):
        self._c = self._c - arg._c
        return self

    def __sub__(self,TailGf y):
        c = self.copy()
        c -= y
        return c

    def __imul__(self,arg):
        """ If arg is a scalar, simple scalar multiplication
            If arg is a GF (any object with _data and _tail as in GF), they it is a matrix multiplication, slice by slice
        """
        n = type(arg).__name__
        if n == 'TailGf' :
            self._c = self._c * (<TailGf?>arg)._c
        elif n in ['float','int', 'complex'] : 
            self._c = as_dcomplex(arg)* self._c
        else : 
            raise RuntimeError, " argument type not recognized in imul for %s"%arg
        return self

    def __mul_impl__(self, arg, s) : 
        cdef TailGf res = self.copy()
        n = type(arg).__name__
        cdef matrix_view [dcomplex,COrder] a 
        if n == 'TailGf' :
            res._c =  self._c * (<TailGf?>arg)._c
        elif n in ['float','int', 'complex'] : 
            res._c = as_dcomplex(arg) * self._c
        else : 
            a= matrix_view[dcomplex,COrder](matrix[dcomplex,COrder](numpy.array(arg, self.dtype)))
            #res._c =  a * self._c  if s else self._c *a
        return res

    def __mul__(self,arg):
        s = type(self).__name__ != 'TailGf' 
        return self.__mul_impl__(arg, s) if not s else arg.__mul_impl__(self, s)

    def __idiv__(self,arg):
        cdef TailGf me = self
        me._c = me._c / as_dcomplex(arg)
        return self

    def __div_impl_(self, arg, s):
        if s : raise RuntimeError, "Can not divide by a TailGf"
        cdef TailGf res = self.copy()
        if type(arg).__name__  in ['float','int', 'complex'] : 
            res._c = self._c / as_dcomplex(arg)
        else : 
            raise RuntimeError, " argument type not recognized for %s"%arg
        return res

    def __div__(self,arg):
        assert type(self).__name__ == 'GfImFreq' 
        s = type(self).__name__ != 'GfImFreq' 
        return self.__div_impl_(arg, s) if not s else arg.__div_impl_(self, s)
    
    #---- other operations ----
    def zero(self) : 
        """Sets the expansion to 0"""
        self.data[:,:,:] =0

    def transpose (self) : 
        """Transpose the array : new view as in numpy"""
        assert(0)

    def conjugate(self) : 
        """Transpose the array : new view as in numpy"""
        assert(0)
        # Hum, a pb here : shall we return a new object (view) or do it one site
        # ? (BAD !!).
        
    def __write_hdf5__ (self, gr , char * key) :
        h5_write (make_h5_group_or_file(gr), key, self._c)

#----------------  Reading from h5 ---------------------------------------

def h5_read_TailGf( gr, std_string key) : 
    return make_TailGf( h5_extractor[tail]()(make_h5_group_or_file(gr),key))

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (TailGf, read_fun = h5_read_TailGf)

#-----------------------------------------------------
# C -> Python 
#-----------------------------------------------------

cdef inline make_TailGf ( tail x) : 
    return TailGf( data = x.data_view().to_python(), order_min = x.order_min() ) #C_Object = encapsulate (&x))



