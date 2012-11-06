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
from ArrayViewWithIndexConverter import ArrayViewWithIndexConverter,_IndicesConverter
import lazy_expressions, Descriptors
from pytriqs.Base.Utility.myUtils import *
from pytriqs.Base.Plot.protocol import clip_array
from TailGF import TailGF 
from types import *
import GF_Initializers
import lazy_expressions,Descriptors

class GfLocalGeneral :

    def _prepare_init(self,d) :
        if 'Indices' in d: 
            indL = list(d.pop('Indices') )
            indR = indL
        elif 'IndicesL' in d and 'IndicesR' in d :
            indL,indR = list(d.pop('IndicesL')) , list(d.pop('IndicesR'))
        else: 
            raise ValueError, "No Indices !!"
        
        self.N1 = len(indL)
        self.N2 = len(indR)
        self._IndicesL = indL
        self._IndicesR = indR
        self._myIndicesGFBlocL = _IndicesConverter(indL)
        self._myIndicesGFBlocR = _IndicesConverter(indR)

        self._tail = d.pop('Tail', None)
        if not self._tail : self._tail = TailGF(OrderMin=-1, size=10, IndicesL=indL, IndicesR=indR)

        self.__data_raw = d.pop('Data', None)
        if self.__data_raw == None: self.__data_raw = numpy.zeros((len(indL),len(indR),len(d['Mesh'])), numpy.complex) #, order='F')

        self._mesh = d.pop('Mesh')
        self.Name = d.pop('Name','g')
        self.Note = d.pop('Note','') 
        
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 
        return  (self._mesh, self.__data_raw, self._tail)

    #-------------   COPY ----------------------------------------

    def copy (self) : 
        new_g = self.__class__(IndicesL = self._IndicesL,
                               IndicesR = self._IndicesR,
                               Mesh = self._mesh,
                               Name = self.Name, Note = self.Note)
        new_g._tail.copyFrom(self._tail)
        new_g._data.array[:,:,:] = self._data.array[:,:,:]
        return new_g
        
    # Should I do more compatibility checks?
    def copyFrom(self, G) :
        self._data.array[:,:,:] = G._data.array[:,:,:]
        self._tail.copyFrom(G._tail)

    #------------------  Access to data ---------------------------------------
    @property
    def mesh(self):
             return self._mesh  

    @property
    def _data(self) : 
        return ArrayViewWithIndexConverter(self.__data_raw, self._IndicesL, self._IndicesR, None)

    @_data.setter
    def _data(self, value) : self.__data_raw[:,:,:] = value 

    @property
    def _data_raw(self) : 
        return self.__data_raw

    @_data.setter
    def _data_raw(self, value) : self.__data_raw[:,:,:] = value 

    #-------------  Indices management ---------------------

    @property
    def Indices(self) : 
        """A generator of the indices"""
        if self._IndicesL != self._IndicesR : raise RuntimeError, "Indices R and L are not the same. I can not give you the Indices"
        for ind in self._IndicesL: 
            yield ind

    @property
    def IndicesL(self) : 
        for ind in self._IndicesL: 
            yield ind
    
    @property
    def IndicesR(self) : 
        for ind in self._IndicesR: 
            yield ind
    
    def arrayWithIndices(self) :
        """ 
        Returns ArrayViewWithIndexConverter with : 
            * the indices of the Green function
            * the correct size
            * an array initialized to 0
        """
        return ArrayViewWithIndexConverter(A = numpy.zeros((self.N1,self.N2), dtype = numpy.complex_),
                                           IndicesL = self._IndicesL, IndicesR = self._IndicesR)

    #---------------------   [  ] operator        ------------------------------------------
    
    def _make_slice(self,sl1, sl2): 
        return self.__class__(IndicesL = self._IndicesL[sl1],
                              IndicesR = self._IndicesR[sl2],
                              Name = self.Name,
                              Mesh = self._mesh,
                              Data = self._data.array[sl1,sl2,:],
                              Tail = self._tail._make_slice(sl1,sl2))

    def __getitem__(self,key):
        """Key is a tuple of index (n1,n2) as defined at construction"""
        if len(key) !=2 : raise ValueError, "[ ] must be given two arguments"
        try :
            sl1 = self._myIndicesGFBlocL.convertToNumpyIndex(key[0])
        except IndexError, ValueError:
            raise IndexError, "Indices %s incorrect. Indices are %s"%(key[0],list(self._IndicesL))

        try :
            sl2 = self._myIndicesGFBlocR.convertToNumpyIndex(key[1])
        except IndexError, ValueError:
            raise IndexError, "Indices %s incorrect. Indices are %s"%(key[1],list(self._IndicesR))

        return self._make_slice (sl1, sl2)

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

    def __reduce__(self):
        return call_factory_from_dict, (self.__class__,self.__reduce_to_dict__())
        
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

    # a classmethod receives the names of the class instead of self.
    @classmethod
    def __factory_from_dict__(CLS,value):
        # a little treatment for backward compatibility
        for s in [ 'Indices' ,'IndicesR' ,'IndicesL' ] : 
            if s in value and  type(value[s]) == type(''):  value[s] = eval(value[s])
        return CLS(**value)
        
    #---------------- Repr, str ---------------------------------

    def __str__ (self) : 
        return self.Name if self.Name else repr(self)
    
    def __repr__(self) : 
        return """%s %s : IndicesL = %s, IndicesR = %s"""%(self.__class__.__name__, self.Name,
          [x for x in self.IndicesL], [x for x in self.IndicesR])

    #--------------   PLOT   ---------------------------------------

    class _Plot_Wrapper_Partial_Reduce : 
        """ Internal Use"""
        def __init__(self, obj,  **opt) : 
            self.obj, self.opt = obj,opt
        def _plot_(self,Options) : 
            Options.update(self.opt)
            return self.obj._plot_(Options)
    
    @property
    def real(self):
        "Use self.real in a plot to plot only the real part"
        return _Plot_Wrapper_Partial_Reduce(self,RI='R')

    @property
    def imag(self): 
        "Use self.imag in a plot to plot only the imag part"
        return _Plot_Wrapper_Partial_Reduce(self,RI='I')
    
    #----------------

    def _plot_base (self, OptionsDict, xlabel, ylabel, use_ris, X):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        Name = OptionsDict.pop('Name', '' )  # consume it
        NamePrefix = OptionsDict.pop('NamePrefix', '' )  # consume it
        if Name and NamePrefix : raise ValueError, 'Name and NamePrefix can not be used at the same time'
        if NamePrefix : name_save, self.Name = self.Name, Name or NamePrefix

        rx = OptionsDict.pop('x_window',None ) # consume it
        sl = clip_array (X, *rx) if rx else slice(len(X)) # the slice due to clip option x_window

        def mdic( prefix, f) : 
           return [{'Type' : "XY", 
                    'xlabel' : xlabel,
                    'ylabel' : ylabel (self.Name),
                    'xdata' : X[sl],
                    'label' : Name if Name else prefix + B.Name ,
                    'ydata' : f( B._data.array[0,0,sl] ) } for (i,j,B) in self ] 
    
        if use_ris : 
            ris = OptionsDict.pop('RI','RI') 
            if   ris == "R" : 
                res = mdic( 'Re ', lambda x : x.real)
            elif ris == "I" : 
                res = mdic( 'Im ', lambda x : x.imag)
            elif ris == "S" :
                res = mdic( '', lambda x : -1/numpy.pi *x.imag)
            elif ris == 'RI' :
                 res = mdic( 'Re ', lambda x : x.real) + mdic( 'Im ', lambda x : x.imag)
            else : 
                 raise ValueError, "RIS flags meaningless %s"%ris
        else: 
            res = mdic( '', lambda x : x)
            
        if NamePrefix: self.Name = name_save
        return res 
 
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
        X = [x.imag for x in self._mesh] if self._mesh.TypeGF == GF_Type.Imaginary_Frequency  else [x for x in self._mesh]
        X, data = numpy.array(X), self._data.array
        if x_window :
          sl = clip_array (X, *x_window) if x_window else slice(len(X)) # the slice due to clip option x_window
          X, data = X[sl],  data[:,:,sl]
        if flatten_y and data.shape[:2]==(1,1) : data = data[0,0,:]
        return X, data
        
    #--------  LAZY expression system -----------------------------------------

    def __lazy_expr_eval_context__(self) : 

        def _is_compatible_for_ops(self, g): 
            m1,m2  = self._mesh, g._mesh
            return m1 is m2 or m1 == m2

        class ctx : 
            def __init__ (self, G) : 
                self.G = G
            def __eq__ (self, y) :
                return isinstance(y, self.__class__) and self.G._is_compatible_for_ops(y.G)
            def __call__ (self, x) : 
                if not isinstance(x, Descriptors.Base) : return x
                tmp = self.G.copy()
                x(tmp)
                return tmp
        
        return ctx(self)

    def __le__(self, other) :
        """ Forbidden : to avoid typo with <<="""
        raise RuntimeError, " Operator <= not defined "
    
    def __ilshift__(self, A): 
        """ A can be two things :
          * G <<= any_GF_Initializers will init the GFBloc with the initializer
          * G <<= g2 where g2 is a GFBloc will copy g2 into self
        """
        if isinstance(A, self.__class__) : 
            if self is not A : self.copyFrom(A) # otherwise it is useless AND does not work !!
        elif isinstance(A, lazy_expressions.lazy_expr) : # A is a lazy_expression made of GF, scalars, descriptors 
            A2= Descriptors.convert_scalar_to_Const(A)
            def e_t (x) : 
                if not isinstance(x, Descriptors.Base) : return x
                tmp = self.copy()
                x(tmp)
                return tmp
            #e_t2 = self.__lazy_expr_eval_context__()
            self.copyFrom ( lazy_expressions.eval_lazy_expr(e_t, A2) )
        elif isinstance(A, lazy_expressions.lazy_expr_terminal) : #e.g. g<<= SemiCircular (...) 
            self <<= lazy_expressions.lazy_expr(A)
        elif Descriptors.is_scalar(A) : #in the case it is a scalar .... 
            self <<= lazy_expressions.lazy_expr(A)
        elif isinstance(A, GF_Initializers.Base) : # backwards compatibility, deprecated
            A(self)
        else :
            raise RuntimeError, " GF Block : <<= operator : RHS not understood"
        return self

    #--------------------  Arithmetic operations  ---------------------------------

    def __iadd__(self,arg):
        d,t = self._data.array, self._tail
        if hasattr(arg,"_data") : # a GF
            d[:,:,:] += arg._data.array
            t += arg._tail
        elif isinstance(arg,numpy.ndarray): # an array considered as a constant function 
            for om in range (d.shape[-1]) : d[:,:,om ] += arg
            t[0].array[:,:] += arg
        elif isinstance(arg,list): # a list
            arg = numpy.array(arg)
            for om in range (d.shape[-1]) : d[:,:,om ] += arg
            t[0].array[:,:] += arg
        elif Descriptors.is_scalar(arg): # just a scalar
            arg = arg*numpy.identity(self.N1)
            for om in range (d.shape[-1]) : d[:,:,om ] += arg
            t[0].array[:,:] += arg
        else:
            raise NotImplementedError
        return self

    def __add__(self,y):
        c = self.copy()
        try: c += y
        except NotImplementedError: return NotImplemented
        return c

    def __radd__(self,y): return self.__add__(y)

    def __isub__(self,arg):
        d,t = self._data.array, self._tail
        if hasattr(arg,"_data") : # a GF
            d[:,:,:] -= arg._data.array
            t -= arg._tail
        elif isinstance(arg,numpy.ndarray): # an array considered as a constant function 
            for om in range (d.shape[-1]) : d[:,:,om ] -= arg
            t[0].array[:,:] -= arg
        elif isinstance(arg,list): # a list
            arg = numpy.array(arg)
            for om in range (d.shape[-1]) : d[:,:,om ] -= arg
            t[0].array[:,:] -= arg
        elif Descriptors.is_scalar(arg): # just a scalar
            arg = arg*numpy.identity(self.N1)
            for om in range (d.shape[-1]) : d[:,:,om ] -= arg
            t[0].array[:,:] -= arg
        else:
          raise NotImplementedError
        return self

    def __sub__(self,y):
        c = self.copy()
        try: c -= y
        except NotImplementedError: return NotImplemented
        return c

    def __rsub__(self,y):
        c = (-1)*self.copy()
        try: c += y
        except NotImplementedError: return NotImplemented
        return c

    def __imul__(self,arg):
        """ If arg is a scalar, simple scalar multiplication
            If arg is a GF (any object with _data and _tail as in GF), they it is a matrix multiplication, slice by slice
        """
        d,t = self._data.array, self._tail
        if hasattr(arg,"_data") : 
            d2 = arg._data.array
            assert d.shape == d2.shape ," Green function block multiplication with arrays of different size !"
            for om in range (d.shape[-1]) : 
                d[:,:,om ] = numpy.dot(d[:,:,om], d2[:,:,om])
            t *= arg._tail
        elif Descriptors.is_scalar(arg): # a scalar
            d[:,:,:] *= arg
            # to be simplified when the *= scalar for tail will be added !
            for n in range(t.OrderMin,t.OrderMax+1):
                t[n].array[:,:] *= arg
        else:
          print type(arg)
          raise NotImplementedError
        return self

    def __mul__(self,y):
        if hasattr(y,"_data") :
            c = self.copy_with_new_stat(GF_Statistic.Boson if self._mesh.Statistic == y._mesh.Statistic else GF_Statistic.Fermion)
        else:
            c = self.copy()
        try: 
            c *= y
        except NotImplementedError: return NotImplemented
        return c

    def __rmul__(self,x): return self.__mul__(x)

    def __idiv__(self,arg):
        d,t = self._data.array, self._tail
        if hasattr(arg,"_data") :
            d2 = arg._data.array
            assert d.shape == d2.shape ," Green function block multiplication with arrays of different size !"
            for om in range (d.shape[-1]) : 
                d[:,:,om ] = numpy.dot(d[:,:,om], numpy.linalg.inv(d2[:,:,om]))
            t /= arg._tail
        elif Descriptors.is_scalar(arg): # a scalar
            d[:,:,:] /= arg
            for n in range(t.OrderMin,t.OrderMax+1):
                t[n].array[:,:] /= arg
            #t = self._tail; t /= arg
        else:
          raise NotImplementedError
        return self

    def __div__(self,y):
        c = self.copy()
        c /= y
        return c

    #--------------   OTHER OPERATIONS -----------------------------------------------------

    def from_L_G_R (self, L,G,R):
        """ For all argument, replace the matrix by L *matrix * R"""
        d,dg = self._data.array,G._data.array
        for om in range (d.shape[-1]) :
            d[:,:,om ] = numpy.dot(numpy.dot(L,dg[:,:,om]), R)
        self._tail.from_L_T_R(L,G._tail,R)

    def invert(self) : 
        """Invert the matrix for all arguments"""
        d = self._data.array
        for om in range (d.shape[-1]) : 
            d[:,:,om ] = numpy.linalg.inv(d[:,:,om])
        self._tail.invert()

    def replaceByTail(self,start) : 
        d = self._data.array
        t = self._tail
        for n, om in enumerate(self._mesh) : # not the most efficient ...
            if n >= start : d[:,:,n] = t(om).array

    def transpose(self):
        """Transposes the GF Bloc: return a new transposed view"""
        
        return self.__class__( 
                Indices = list(self.Indices),
                Mesh  = self._mesh,
                Data = self._data.array.transpose( (1,0,2) ), 
                Tail = self._tail.transpose(),
                Name = self.Name+'(t)', 
                Note = self.Note)

    def conjugate(self):
        """Complex conjugate of the GF Bloc. It follow the policy of numpy and
        make a copy only if the Green function is complex valued"""
        
        return self.__class__( 
                Indices = list(self.Indices),
                Mesh  = self._mesh,
                Data = self._data.conjugate(), 
                Tail = self._tail.conjugate(self._mesh.TypeGF==GF_Type.Imaginary_Frequency),
                Name = self.Name+'*', 
                Note = self.Note)
    
    #------------------  Density -----------------------------------
    
    def total_density(self) :
        """Trace density"""
        return numpy.trace(self.density())

