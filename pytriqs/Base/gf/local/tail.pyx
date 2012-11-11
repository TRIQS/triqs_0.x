cdef class TailGf:
    """ 
        Doc to be written 
    """
    
    cdef tail _c
    
    def __init__(self, **d):
        """
        """
        self._indL, self._indR = d['IndicesL'], d['IndicesR']
        N1,N2 = len(self._indL), len(self._indR)
        a = d['array'] if 'array' in d else numpy.zeros((N1,N2,d['size'])  ,numpy.complex) #,order='F')
        omin = d['OrderMin']
        self._c =  tail( array_view[dcomplex,THREE,COrder](a), omin)
    
    def invert(self) :
        self._c = 1.0/self._c
    
    def __fill_a( self,a, value) : 
        if hasattr(value,'shape') : 
            if a.shape[:2] != value.shape[:2] : 
                raise RuntimeError, "shape mismatch"
            m = min(value.shape[2],a.shape[2]) 
            if a.shape[2] > value.shape[2] : 
                a[:,:,m:] =0
            a[:,:, 0:m] = value[:,:, 0:m]
        else : 
            a = value

    property OrderMin : 
        """Min order of the expansion"""
        def __get__(self) : return self._c.order_min()

    property OrderMax : 
        """Max order of the expansion"""
        def __get__(self) : return self._c.order_max()

    property size : 
        """Length of the expansion"""
        def __get__(self) : return self._c.size()

    property _data_raw : 
        """Access to the data array"""
    
        def __get__(self) : 
            return self._c.data_view().to_python()

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            self.__fill_a(a,value)   
    
    property _data : 
        """Access to the data array"""
    
        def __get__(self) : 
            return ArrayViewWithIndexConverter(self._c.data_view().to_python(), self._indL, self._indR, None)

        def __set__ (self, value) :
            cdef object a = self._c.data_view().to_python()
            self.__fill_a(a,value)   

    def _make_slice(self, sl1, sl2):
        return self.__class__(IndicesL = self._indL[sl1], IndicesR = self._indR[sl2], array = self._data_raw[sl1,sl2,:], OrderMin = self._omin)
        #return self.__class__(IndicesL = self._indL[sl1], IndicesR = self._indR[sl2], array = numpy.asfortranarray(self._data_raw[sl1,sl2,:]), OrderMin = self._omin)

    def copy(self) : 
        #return self.__class__(IndicesL = self._indL, IndicesR = self._indR, array = self._data_raw.copy(order='F'), OrderMin = self._omin)
        return self.__class__(IndicesL = self._indL, IndicesR = self._indR, array = self._data_raw.copy(order='C'), OrderMin = self._omin)

    # Should I do more compatibility checks?
    def copyFrom(self, T) : 
        self._omin = T._omin
        self._data_raw = T._data_raw.copy(order='C')

    def __repr__ (self) :
        return string.join([ "%s"%self[r].array + (" /" if r<0 else "") + " Om^%s"%(abs(r)) for r in range(self.OrderMin, self.OrderMax+1) ] , " + ")

    def __getitem__(self,i) :
        """Returns the i-th coefficient of the expansion, or order Om^i"""
        if not self.has_coef(i) : raise IndexError, "Index %s is out of range"%i
        return ArrayViewWithIndexConverter(self._data_raw[:,:,i-self.OrderMin], self._indL, self._indR)

    def __setitem__(self,i, val) :
        """Sets the i-th coefficient of the expansion, or order Om^i"""
        if not self.has_coef(i) : raise IndexError, "Index %s is out of range"%i
        self._data_raw[:,:,i-self._omin] = val

    def has_coef(self, i):
        return (i >= self.OrderMin) and (i <= self.OrderMax)

    def __call__(self, n) :
        if not self.has_coef(n) : raise IndexError, "Index %s is out of range"%n
        return self._c(n).to_python()   

    #-----------------------------------------------------
    def __reduce__(self):
        return call_factory_from_dict, (self.__class__,self.__reduce_to_dict__())
  
    def __reduce_to_dict__(self):
        indL = repr(tuple(self._indL))
        indR = repr(tuple(self._indR))
        assert(eval(indL)==tuple(self._indL))
        assert(eval(indR)==tuple(self._indR))
        return {'IndicesL' : indL,
                'IndicesR' : indR,
                'array' : self._data_raw,
                'OrderMin' : self.OrderMin,
                'size' : self.size }

    # a classmethod receives the names of the class instead of self.
    @classmethod
    def __factory_from_dict__(CLS,d):
        # a little treatment for backward compatibility
        for s in [ 'Indices' ,'IndicesR' ,'IndicesL' ] : 
            if s in d and  type(d[s]) == type(''):  d[s] = eval(d[s])
        return CLS(**d)
 
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
            res._c =  a * self._c  if s else self._c *a
        return res

    def __mul__(self,arg):
        s = type(self).__name__ != 'TailGf' 
        return self.__mul_impl__(self, arg, s) if not s else self.__mul_impl__(self, arg, s)

    def __idiv__(self,arg):
        cdef TailGf me = self
        me._c = me._c / as_dcomplex(arg)
        return self

    def __div_impl_(self, arg, s):
        if s : raise RuntimeError, "Can not divide by an TailGf"
        cdef TailGf res = self.copy()
        if type(arg).__name__  in ['float','int', 'complex'] : 
            res._c = self._c / as_dcomplex(arg)
        else : 
            raise RuntimeError, " argument type not recognized for %s"%arg
        return res

    def __div__(self,arg):
        s = type(self).__name__ != 'GfImFreq' 
        return self.__div_impl__(self, arg, s) if not s else self.__div_impl__(self, arg, s)

    #---- other operations ----
    def zero(self) : 
        """Sets the expansion to 0"""
        self._data_raw[:,:,:] =0

    def transpose (self) : 
        """Transpose the array : new view as in numpy"""
        assert(0)

    def conjugate(self) : 
        """Transpose the array : new view as in numpy"""
        assert(0)
        # Hum, a pb here : shall we return a new object (view) or do it one site
        # ? (BAD !!).

#-----------------------------------------------------
# C -> Python 
#-----------------------------------------------------

cdef inline make_TailGf ( tail x) : 
    return TailGf(C_Object = encapsulate (&x))

#-----------------------------------------------------
#  Register the class for HDF_Archive
#-----------------------------------------------------

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (TailGf)



