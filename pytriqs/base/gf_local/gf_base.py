
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
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

from pytriqs_GF import GF_Statistic,GF_Type,TailGf,MeshGf
from pytriqs.base.utility.my_utils import *
import numpy
from types import *
from array_view import ArrayViewWithIndexConverter,_IndicesConverter
import gf_init
import lazy_expressions, descriptors
from pytriqs.base.plot.protocol import clip_array

class _Plot_Wrapper_Partial_Reduce : 
    """ Internal Use"""
    def __init__(self, obj,  **opt) : 
        self.obj, self.opt = obj,opt
    def _plot_(self, options) : 
        options.update(self.opt)
        return self.obj._plot_(options)

class GfBase:
    
    def _init_base__(self,d) :
        if 'indices' in d: 
            indL = list(d.pop('indices') )
            indR = indL
        elif 'indicesL' in d and 'indicesR' in d :
            indL,indR = list(d.pop('indicesL')) , list(d.pop('indicesR'))
        else: 
            raise ValueError, "No Indices !!"
        
        self._myIndicesGFBlocL = _IndicesConverter(indL)
        self._myIndicesGFBlocR = _IndicesConverter(indR)
        mesh = d.pop('mesh')
        data = d.pop('data', None)
        tail = d.pop('tail', None)
        if not tail : tail = TailGf(-2,10,indL,indR)
        self.name = d.pop('name','g')
        self.note = d.pop('note','') 
        
        self._param_for_cons = (indL,indR,data,mesh,tail)
        assert len(d) ==0, "Unknown parameters in BlockGf constructions %s"%d.keys() 

    #---------------------------------------------------------------------------------
 
    @property
    def _data(self) : 
        return ArrayViewWithIndexConverter(self._data_c_array, self._IndicesL, self._IndicesR, None)

    @_data.setter
    def _data(self, value) : self._data_c_array[:,:,:] = value 
    
    #---------------------------------------------------------------------------------
    
    def _make_slice(self,sl1, sl2): 
        return self.__class__(indicesL = self._IndicesL[sl1],
                              indicesR = self._IndicesR[sl2],
                              name = self.name,
                              mesh = self.mesh,
                              data = self._data.array[sl1,sl2,:],
                              tail = TailGf(self._tail,sl1,sl2))

    #-----------------------------------------------------

    def copy (self) : 
        new_g = self.__class__(indicesL = self._IndicesL,
                               indicesR = self._IndicesR,
                               mesh = self.mesh,
                               name = self.name, note = self.note)
        new_g._tail.copy_from(self._tail)
        new_g._data.array[:,:,:] = self._data.array[:,:,:]
        return new_g
        
    #-----------------------------------------------------
        
    def copy_with_new_stat(self,stat) :
        new_g = self.__class__(indicesL = self._IndicesL,
                               indicesR = self._IndicesR,
                               mesh = MeshGf(self.mesh,stat),
                               name = self.name, note = self.note)
        new_g._tail.copy_from(self._tail)
        new_g._data.array[:,:,:] = self._data.array[:,:,:]
        return new_g        

    #-----------------------------------------------------

    def _is_compatible_for_ops(self, g): 
        m1,m2  = self.mesh, g.mesh
        return m1 is m2 or m1 == m2

    #-----------------------------------------------------

    def __reduce_to_dict__(self):
        indL = repr(tuple(self.indicesL))
        indR = repr(tuple(self.indicesR))
        assert(eval(indL)==tuple(self.indicesL))
        assert(eval(indR)==tuple(self.indicesR))
        return {'IndicesL' : indL,
                'IndicesR' : indR,
                'Data' : self._data.array,
                'Mesh' : self.mesh,
                'Tail' : self._tail,
                'Name' : self.name,
                'Note' : self.note }

    # a classmethod receives the names of the class instead of self.
    @classmethod
    def __factory_from_dict__(CLS,value):
        # a little treatment for backward compatibility
        for s in [ 'Indices' ,'IndicesR' ,'IndicesL' ] : 
            if s in value and  type(value[s]) == type(''):  value[s] = eval(value[s])
        return CLS(indicesL = value['IndicesL'],
                   indicesR = value['IndicesR'],
                   data = value['Data'],
                   mesh = value['Mesh'],
                   tail = value['Tail'],
                   name = value['Name'],
                   note = value['Note'])
 
    #-------------------------------------------------
                
    def __lazy_expr_eval_context__(self) : 

        class ctx : 
            def __init__ (self, G) : 
                self.G = G
            def __eq__ (self, y) :
                return isinstance(y, self.__class__) and self.G._is_compatible_for_ops(y.G)
            def __call__ (self, x) : 
                if not isinstance(x, descriptors.Base) : return x
                tmp = self.G.copy()
                x(tmp)
                return tmp
        
        return ctx(self)

    def __ilshift__(self, A): 
        """ A can be two things :
          * G <<= any_gf_init will init the BlockGf with the initializer
          * G <<= g2 where g2 is a BlockGf will copy g2 into self
        """
        if isinstance(A, self.__class__) : 
            if self is not A : self.copy_from(A) # otherwise it is useless AND does not work !!
        elif isinstance(A, lazy_expressions.LazyExpr) : # A is a lazy_expression made of BlockGf, scalars, descriptors 
            A2= descriptors.convert_scalar_to_const(A)
            def e_t (x) : 
                if not isinstance(x, descriptors.Base) : return x
                tmp = self.copy()
                x(tmp)
                return tmp
            #e_t2 = self.__lazy_expr_eval_context__()
            self.copy_from ( lazy_expressions.eval_lazy_expr(e_t, A2) )
        elif isinstance(A, lazy_expressions.LazyExprTerminal) : #e.g. g<<= SemiCircular (...) 
            self <<= lazy_expressions.LazyExpr(A)
        elif descriptors.is_scalar(A) : #in the case it is a scalar .... 
            self <<= lazy_expressions.LazyExpr(A)
        elif isinstance(A, gf_init.Base) : # backwards compatibility, deprecated
            A(self)
        else :
            raise RuntimeError, " BlockGf: <<= operator : RHS not understood"
        return self

    #-----------------------------------------------------
    
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
        X = [x.imag for x in self.mesh] if self.mesh.TypeGF == GF_Type.Imaginary_Frequency  else [x for x in self.mesh]
        X, data = numpy.array(X), self._data.array
        if x_window :
          sl = clip_array (X, *x_window) if x_window else slice(len(X)) # the slice due to clip option x_window
          X, data = X[sl],  data[:,:,sl]
        if flatten_y and data.shape[:2]==(1,1) : data = data[0,0,:]
        return X, data

    #--------------------  Arithmetic operations  ---------------------------------

    def __iadd__(self,arg):
        d,t = self._data.array, self._tail
        if hasattr(arg,"_data") : # a Green's function
            d[:,:,:] += arg._data.array
            t += arg._tail
        elif isinstance(arg,numpy.ndarray): # an array considered as a constant function 
            for om in range (d.shape[-1]) : d[:,:,om ] += arg
            t[0].array[:,:] += arg
        elif isinstance(arg,list): # a list
            arg = numpy.array(arg)
            for om in range (d.shape[-1]) : d[:,:,om ] += arg
            t[0].array[:,:] += arg
        elif descriptors.is_scalar(arg): # just a scalar
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
        if hasattr(arg,"_data") : # a Green's function
            d[:,:,:] -= arg._data.array
            t -= arg._tail
        elif isinstance(arg,numpy.ndarray): # an array considered as a constant function 
            for om in range (d.shape[-1]) : d[:,:,om ] -= arg
            t[0].array[:,:] -= arg
        elif isinstance(arg,list): # a list
            arg = numpy.array(arg)
            for om in range (d.shape[-1]) : d[:,:,om ] -= arg
            t[0].array[:,:] -= arg
        elif descriptors.is_scalar(arg): # just a scalar
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
            If arg is a BlockGf(any object with _data and _tail as in BlockGf), they it is a matrix multiplication, slice by slice
        """
        d,t = self._data.array, self._tail
        if hasattr(arg,"_data") : 
            d2 = arg._data.array
            assert d.shape == d2.shape ," Green function block multiplication with arrays of different size !"
            for om in range (d.shape[-1]) : 
                d[:,:,om ] = numpy.dot(d[:,:,om], d2[:,:,om])
            t *= arg._tail
        elif descriptors.is_scalar(arg): # a scalar
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
            c = self.copy_with_new_stat(GF_Statistic.Boson if self.mesh.statistic == y.mesh.statistic else GF_Statistic.Fermion)
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
        elif descriptors.is_scalar(arg): # a scalar
            d[:,:,:] /= arg
            for n in range(t.OrderMax):
                t[n].array[:,:] /= arg
            #t = self._tail; t /= arg
        else:
          raise NotImplementedError
        return self

    def __div__(self,y):
        c = self.copy()
        c /= y
        return c

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

    def replace_by_tail(self,start) : 
        d = self._data.array
        t = self._tail
        for n, om in enumerate(self.mesh) : # not the most efficient ...
            if n >= start : d[:,:,n] = t(om).array

    #-------------------------------------------------------------------

    # remove this
    def copy_and_transpose(self):
        """deprecated :Transposes the BlocGf, but gives back a new instance. Does not change the actual BlockGf!"""
        assert 0
        gtmp = self.copy()
         
        #data:
        M = [x.imag for x in self.mesh]
        for iom in range(len(M)):
            gtmp._data.array[:,:,iom] = self._data.array[:,:,iom].transpose()

        #tails:
        for n in range(self._tail.OrderMax):
            gtmp._tail[n].array[:,:] = self._tail[n].array.transpose()

        return gtmp

    def transpose(self):
        """Transposes the BlocGf: return a new transposed view"""
        
        return self.__class__( 
                indices = list(self.indices),
                mesh  = self.mesh,
                data = self._data.array.transpose( (1,0,2) ), 
                tail = self._tail.transpose(),
                name = self.name+'(t)', 
                note = self.note)

    def conjugate(self):
        """Complex conjugate of the BlocGf. It follow the policy of numpy and
        make a copy only if the Green function is complex valued"""
        
        return self.__class__( 
                indices = list(self.indices),
                mesh  = self.mesh,
                data = self._data.array.conjugate(), 
                tail = self._tail.conjugate(self.mesh.TypeGF==GF_Type.Imaginary_Frequency),
                name = self.name+'*', 
                note = self.note)

    #-----------------------------------------------------

    # Put it out a a free funciton in a module. Nothing to do here...
    # why is this here ???
    def delta(self) :
        """Computes delta from self ...."""
        if self.mesh.TypeGF not in [GF_Type.Real_Frequency, GF_Type.Imaginary_Frequency] :
            raise RuntimeError, "Can not compute delta for this BlockGf"
        G0 = self if self._tail.OrderMin <=-1 else inverse(self)
        tmp = G0.copy()
        tmp <<= gf_init.A_Omega_Plus_B(G0._tail[-1], G0._tail[0])
        tmp -= G0
        return tmp

    #-------------  Indices management ---------------------

    def __get_Indices(self) : 
	"""A generator of the indices"""
        if self._IndicesL != self._IndicesR : raise RuntimeError, "Indices R and L are not the same. I can not give you the Indices"
        for ind in self._IndicesL: 
            yield ind

    def __get_IndicesL(self) : 
        for ind in self._IndicesL: 
            yield ind
    
    def __get_IndicesR(self) : 
        for ind in self._IndicesR: 
            yield ind
    
    indices = property(__get_Indices)
    indicesL = property(__get_IndicesL)
    indicesR = property(__get_IndicesR)
            
    def array_with_indices(self) :
        """ 
        Returns ArrayViewWithIndexConverter with : 
            * the indices of the Green function
            * the correct size
            * an array initialized to 0
        """
        return ArrayViewWithIndexConverter(A = numpy.zeros((self.N1,self.N2), dtype = numpy.complex_),
                                           IndicesL = self._IndicesL, IndicesR = self._IndicesR)

    #---------------------   [  ] operator        ------------------------------------------
    
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
 
     #-----------------------------------------------------

    def __reduce__(self):
        return call_factory_from_dict, (self.__class__,self.__reduce_to_dict__())
        
    #-------------------------------------------------

    def __le__(self, other) :
        """ Forbidden : to avoid typo with <<="""
        raise RuntimeError, " Operator <= not defined "
        
    #-------------------------------------------------

    def __str__ (self) : 
	return self.name if self.name else repr(self)
    
    #-------------------------------------------------

    def __repr__(self) : 
	return """%s %s :  Beta = %.3f; IndicesL = %s, IndicesR = %s """%(self.__class__.__name__, self.name,
          self.beta, [x for x in self.indicesL], [x for x in self.indicesR])

    #-------------------------------------------------

    def __iter__(self) :
	for i in self._IndicesL : 
	    for j in self._IndicesR :
                b =self[i,j]
		b.name = "%s_%s_%s"%(self.name if hasattr(self,'name') else '',i,j)
		yield i,j,b

    #-----------------------------------------------------
 
    def __get_real(self) : return _Plot_Wrapper_Partial_Reduce(self,RI='R')
    def __get_imag(self) : return _Plot_Wrapper_Partial_Reduce(self,RI='I')
    real = property(__get_real,None, doc = "Use self.real in a plot to plot only the real part")
    imag = property(__get_imag,None, doc = "Use self.imag in a plot to plot only the imag part")

    #-----------------------------------------------------
    def total_density(self) :
        """Trace density"""
        return numpy.trace(self.density())

    #-----------------------------------------------------

    def _plot_base (self, opt_dict, xlabel, ylabel, use_ris, X):
        """ Plot protocol. opt_dict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        Name = opt_dict.pop('name', '' )  # consume it
        NamePrefix = opt_dict.pop('name_prefix', '' )  # consume it
        if Name and NamePrefix : raise ValueError, 'name and name_prefix can not be used at the same time'
        if NamePrefix : name_save, self.name = self.name, Name or NamePrefix

        rx = opt_dict.pop('x_window',None ) # consume it
        sl = clip_array (X, *rx) if rx else slice(len(X)) # the slice due to clip option x_window

        def mdic( prefix, f) : 
           return [{'type' : "XY", 
                    'xlabel' : xlabel,
                    'ylabel' : ylabel (self.name),
                    'xdata' : X[sl],
                    'label' : Name if Name else prefix + B.name ,
                    'ydata' : f( B._data.array[0,0,sl] ) } for (i,j,B) in self ] 
    
        if use_ris : 
            ris = opt_dict.pop('RI','RI') 
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
            
        if NamePrefix: self.name = name_save
        return res 
 
#-----------------------------------------------------

from pytriqs.base.archive.hdf_archive_schemes import register_class
register_class (TailGf)
register_class (MeshGf)
