################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011-2012 by M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
import numpy
import lazy_expressions, Descriptors
import pytriqs.Base.Utility.myUtils 
from pytriqs.Base.Plot.protocol import clip_array
import lazy_expressions,Descriptors
from types import IntType,SliceType,StringType
from tools import PlotWrapperPartialReduce, lazy_ctx, IndicesConverter,get_indices_in_dict, py_deserialize
import impl_plot

# Function that transcribe the indices to C++
cdef indices_2_t make_c_indices(indicesL, indicesR) : 
    cdef vector[vector[std_string]] res
    cdef vector[std_string] vl,vr
    for i in indicesL : vl.push_back(i)
    for i in indicesR : vr.push_back(i)
    res.push_back(vl); res.push_back(vr)
    return indices_2_t(res)

cdef class _ImplGfLocal :
    
    cdef object _name, dtype, _mesh, _data, _singularity, _symmetry, _indices
    cdef readonly  _derived

    def __init__(self, mesh, data, singularity, symmetry, indices, name, derived ) : 
        self._mesh, self._data, self._singularity, self._symmetry, self._indices= mesh, data, singularity, symmetry, indices
        self._name = name
        self._derived = derived

    property mesh : 
        """Mesh"""
        def __get__(self): return self._mesh
    
    property tail : 
        def __get__(self): return self._singularity
        def __set__(self,TailGf t): 
            assert (self.N1, self.N2, self._singularity.size) == (t.N1, t.N2, t.size)
            self._singularity.copyFrom (t)

    property data : 
        """Access to the data array"""
        def __get__(self) : return self._data
        def __set__ (self, value) : self._data[:,:,:] = value
   
    property N1 : 
        def __get__(self): return self.data.shape[0]

    property N2 : 
        def __get__(self): return self.data.shape[1]

    property indicesL : 
        """Indices ..."""
        def __get__(self) :
            v = self._indices[0]
            for ind in v:
                yield ind
    
    property indicesR : 
        """Indices ..."""
        def __get__(self) : 
            v = self._indices[1]
            for ind in v:
                yield ind

    property indices : 
        """Indices ..."""
        def __get__(self) : 
            inds =  self._indices
            if inds[0] != inds[1]: raise RuntimeError, "Indices R and L are not the same. I can not give you the Indices"
            for i in inds[0]:
                yield i

    property Name : 
        """Name of the Green function (for plots, etc...) """
        def __get__(self) : return self._name
        def __set__(self,val) : self._name = str(val)

    #-------------- Reduction -------------------------------
    # Incorrect .... ** missing
    def __reduce__(self):
        return self._derived, { 'Mesh' : self._mesh, 'Data' : self._data, 
                'Tail' : self._singularity, 'Symmetry' : self._symmetry,
                'Indices' : self._indices, 'Name' : self._name } 

    #-------------   COPY ----------------------------------------

    def copy (self) : 
        cls, args = self.__reduce__()
        print cls, args
        return cls(**args)
        
    def copy_from(self, X) :
        assert self._derived is X._derived
        assert self._mesh == X.mesh
        self.data = X.data
        self.tail = X.tail
        #assert list(self._indices)== list( X._indices)
        #self._symmetry = X._symmetry
        self._name = X.Name

    #---------------------   [  ] operator        ------------------------------------------
    
    def __getitem__(self,key):
        """Key is a tuple of index (n1,n2) as defined at construction"""
        if len(key) !=2 : raise KeyError, "[ ] must be given two arguments"
        sl1, sl2 = key
        if type(sl1) == StringType and type(sl2) == StringType :
            # Convert the indices to integer
            indices_converter = [ IndicesConverter(self.indicesL), IndicesConverter(self.indicesR)]
            sl1, sl2 =  [ indices_converter[i].convertToNumpyIndex(k) for i,k in enumerate(key) ]
        if type (sl1) != slice : sl1 = slice (sl1, sl1+1)
        if type (sl2) != slice : sl2 = slice (sl2, sl2+1)
        return self.__class__(IndicesL = list(self.indicesL)[sl1],
                              IndicesR = list(self.indicesR)[sl2],
                              Name = self.Name,
                              Mesh = self.mesh,
                              Data = self.data[sl1,sl2,:],
                              Tail = self.tail._make_slice(sl1,sl2))

    def __setitem__(self,key,val):
        g = self.__getitem__(key)
        g <<= val

    #------------- Iteration ------------------------------------

    def __iter__(self) :
        for i in self.indicesL : 
            for j in self.indicesR :
                b =self[i,j]
                b.Name = "%s_%s_%s"%(self.Name if hasattr(self,'Name') else '',i,j)
                yield i,j,b
 
    #---------------- Repr, str ---------------------------------

    def __str__ (self) : 
        return self.Name if self.Name else repr(self)
    
    def __repr__(self) : 
        return """%s %s : IndicesL = %s, IndicesR = %s"""%(self.__class__.__name__, self.Name,
          [x for x in self.indicesL], [x for x in self.indicesR])

    #--------------   PLOT   ---------------------------------------

    property real : 
        """Use self.real in a plot to plot only the real part"""
        def __get__ (self): return PlotWrapperPartialReduce(self,RI='R')

    property imag : 
        """Use self.imag in a plot to plot only the imag part"""
        def __get__ (self): return PlotWrapperPartialReduce(self,RI='I')
  
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
        X = [x.imag for x in self.mesh] if type(self.mesh) == MeshImFreq else [x for x in self.mesh]
        X, data = numpy.array(X), self.data
        if x_window :
          sl = clip_array (X, *x_window) if x_window else slice(len(X)) # the slice due to clip option x_window
          X, data = X[sl],  data[:,:,sl]
        if flatten_y and data.shape[:2]==(1,1) : data = data[0,0,:]
        return X, data
        
    #--------  LAZY expression system -----------------------------------------

    def __lazy_expr_eval_context__(self) : 
        return lazy_ctx(self)

    def __richcmp__(self, other, op) :
        # Shall I implement == ?
        raise RuntimeError, " Operator not defined "
    
    def __ilshift__(self, A): 
        """ A can be two things :
          * G <<= any_GF_Initializers will init the GFBloc with the initializer
          * G <<= g2 where g2 is a GFBloc will copy g2 into self
        """
        if isinstance(A, self.__class__) : 
            if self is not A : self.copy_from(A) # otherwise it is useless AND does not work !!
        elif isinstance(A, lazy_expressions.lazy_expr) : # A is a lazy_expression made of GF, scalars, descriptors 
            A2= Descriptors.convert_scalar_to_Const(A)
            def e_t (x) : 
                if not isinstance(x, Descriptors.Base) : return x
                tmp = self.copy()
                x(tmp)
                return tmp
            self.copy_from ( lazy_expressions.eval_lazy_expr(e_t, A2) )
        elif isinstance(A, lazy_expressions.lazy_expr_terminal) : #e.g. g<<= SemiCircular (...) 
            self <<= lazy_expressions.lazy_expr(A)
        elif Descriptors.is_scalar(A) : #in the case it is a scalar .... 
            self <<= lazy_expressions.lazy_expr(A)
        else :
            raise RuntimeError, " <<= operator : RHS  not understood"
        return self

    #--------------------  Arithmetic operations  ---------------------------------

    def __iadd__(self,arg):
        assert type(arg).__name__ == self._typename, "Can not add a %s to a %s"%(type(self).__name__ , type(arg).__name__ )
        self.data[:,:,:] += arg.data
        self.tail += arg.tail
        return self

    def __isub__(self,arg):
        assert type(arg).__name__ == self._typename, "Can not substract a %s from a %s"%(type(arg).__name__, type(self).__name__) 
        self.data[:,:,:] -= arg.data
        self.tail -= arg.tail
        return self

    def __add__(self,y):
        c = self.copy()
        c += y
        return c

    def __sub__(self,y):
        c = self.copy()
        c -= y
        return c

    def __imul__(self,arg):
        """ If arg is a scalar, simple scalar multiplication
            If arg is a GF (any object with data and tail as in GF), they it is a matrix multiplication, slice by slice
        """
        n = type(arg).__name__
        if n == self._typename:
            d,d2 = self.data, arg.data
            assert d.shape == d2.shape ," Green function block multiplication with arrays of different size !"
            for om in range (d.shape[-1]) : 
                d[:,:,om ] = numpy.dot(d[:,:,om], d2[:,:,om])
            self.tail = arg.tail * self.tail
        elif n in ['float','int', 'complex'] : 
            self.data *= arg 
            self.tail *= arg  
        else : 
            raise RuntimeError, " argument type not recognized in imul for %s"%arg
        return self

    def __mul__(self,arg):
        cdef int i
        a,b = (self, arg) if type(self).__name__ in  ['float','int', 'complex'] else (arg,self)
        res = b.copy()
        res *=a
        return res

    def imatmul_L(self,L):
        cdef int i
        dot = numpy.dot
        assert type(L).__name__ in  ['ndarray', 'matrix'] 
        A = self.data
        for i in range( A.shape[-1]) :
            A[:,:,i] = dot(L,A[:,:,i])
        return self

    def imatmul_R(self,R):
        cdef int i
        dot = numpy.dot
        assert type(R).__name__ in  ['ndarray', 'matrix'] 
        A = self.data
        for i in range( A.shape[-1]) :
            A[:,:,i] = dot(A[:,:,i],R)
        return self

    def imatmul_LR(self,L,R): 
        cdef int i
        dot = numpy.dot
        assert type(R).__name__ in  ['ndarray', 'matrix'] 
        assert type(L).__name__ in  ['ndarray', 'matrix'] 
        A = self.data
        for i in range( A.shape[-1]) :
            A[:,:,i] = dot(L, dot(A[:,:,i],R))
        return self

    # RENAME THIS !
    def from_L_G_R (self, L,G,R):
        """ For all argument, replace the matrix by L *matrix * R"""
        d,dg = self.data,G.data
        for om in range (d.shape[-1]) :
            d[:,:,om ] = numpy.dot(numpy.dot(L,dg[:,:,om]), R)
        self.tail.data = numpy.dot(numpy.dot(L,G.tail), R)
    
    def __idiv__(self,arg):
        """ If arg is a scalar, simple scalar multiplication
        """
        n = type(arg).__name__
        if n in ['float','int', 'complex'] : 
            self.data /= arg 
            self.tail /= arg  
        else : 
            raise RuntimeError, " argument type not recognized in imul for %s"%arg
        return self

    def __div__(self,arg):
        assert type(arg).__name__ in  ['float','int', 'complex'], "Error in /"
        res = self.copy()
        res /= arg
        return res
        
    # VERY BAD : what if I make a product in omega ? Legendre ?
    #def __mul__(self,y):
    #    if hasattr(y,"_data") :
    #        c = self.copy_with_new_stat(GF_Statistic.Boson if
    #        self._mesh.Statistic == y._mesh.Statistic else
    #        GF_Statistic.Fermion#)
    #    else:
    #        c = self.copy()
    #    try: 
    #        c *= y
    #    except NotImplementedError: return NotImplemented
    #    return c

    #---------------------------------------------------

    def invert(self) : 
        """Invert the matrix for all arguments"""
        pass
        #self._c = inverse_c (self._c)

    #---------------------------------------------------    
    def transpose(self):
        """Transposes the GF Bloc: return a new transposed view"""
       
       ### WARNING : this depends on the C++ layering ....
        return self.__class__( 
                Indices = list(self.indices),
                Mesh  = self.mesh,
                Data = self.data.transpose( (1,0,2) ), 
                Tail = self.tail.transpose(),
                Name = self.Name+'(t)', 
                Note = self.Note)

    #---------------------------------------------------

    def conjugate(self):
        """Complex conjugate of the GF Bloc. It follow the policy of numpy and
        make a copy only if the Green function is complex valued"""
        
        return self.__class__( 
                Indices = list(self.indices),
                Mesh  = self.mesh,
                Data = self.data.conjugate(), 
                Tail = self.tail.conjugate(), 
                Name = self.Name+'*', 
                Note = self.Note)
    
    #------------------  Density -----------------------------------
    
    def total_density(self) :
        """Trace density"""
        return numpy.trace(self.density())

