
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

from pytriqs_GF import GF_Statistic,GF_Type,TailGF,MeshGF
from pytriqs.base.Utility.myUtils import *
import numpy
from types import *
from ArrayViewWithIndexConverter import ArrayViewWithIndexConverter,_IndicesConverter
import GF_Initializers
import lazy_expressions,Descriptors
from pytriqs.base.Plot.protocol import clip_array

class _GFBloc_base_data_tail  :
    
    def _init_base__(self,d) :
        if 'Indices' in d: 
            indL = list(d.pop('Indices') )
            indR = indL
        elif 'IndicesL' in d and 'IndicesR' in d :
            indL,indR = list(d.pop('IndicesL')) , list(d.pop('IndicesR'))
        else: 
            raise ValueError, "No Indices !!"
        
        self._myIndicesGFBlocL = _IndicesConverter(indL)
        self._myIndicesGFBlocR = _IndicesConverter(indR)
        Mesh = d.pop('Mesh')
        Data = d.pop('Data', None)
        Tail = d.pop('Tail', None)
        if not Tail : Tail = TailGF(-2,10,indL,indR)
        self.Name = d.pop('Name','g')
        self.Note = d.pop('Note','') 
        
        self._param_for_cons = (indL,indR, Data,Mesh,Tail)
        assert len(d) ==0, "Unknown parameters in GFBloc constructions %s"%d.keys() 

    #---------------------------------------------------------------------------------
 
    @property
    def _data(self) : 
        return ArrayViewWithIndexConverter(self._data_c_array, self._IndicesL, self._IndicesR, None)

    @_data.setter
    def _data(self, value) : self._data_c_array[:,:,:] = value 
    
    #---------------------------------------------------------------------------------
    
    def _make_slice(self,sl1, sl2): 
        return self.__class__(IndicesL = self._IndicesL[sl1],
                              IndicesR = self._IndicesR[sl2],
                              Name = self.Name,
                              Mesh = self.mesh,
                              Data = self._data.array[sl1,sl2,:],
                              Tail = TailGF(self._tail,sl1,sl2))

    #-----------------------------------------------------

    def copy (self) : 
        new_g = self.__class__(IndicesL = self._IndicesL,
                               IndicesR = self._IndicesR,
                               Mesh = self.mesh,
                               Name = self.Name, Note = self.Note)
        new_g._tail.copyFrom(self._tail)
        new_g._data.array[:,:,:] = self._data.array[:,:,:]
        return new_g
        
    #-----------------------------------------------------
        
    def copy_with_new_stat(self,stat) :
        new_g = self.__class__(IndicesL = self._IndicesL,
                               IndicesR = self._IndicesR,
                               Mesh = MeshGF(self.mesh,stat),
                               Name = self.Name, Note = self.Note)
        new_g._tail.copyFrom(self._tail)
        new_g._data.array[:,:,:] = self._data.array[:,:,:]
        return new_g        

    #-----------------------------------------------------

    def _is_compatible_for_ops(self, g): 
        m1,m2  = self.mesh, g.mesh
        return m1 is m2 or m1 == m2

    #-----------------------------------------------------

    def __reduce_to_dict__(self):
        indL = repr(tuple(self.IndicesL))
        indR = repr(tuple(self.IndicesR))
        assert(eval(indL)==tuple(self.IndicesL))
        assert(eval(indR)==tuple(self.IndicesR))
        return {'IndicesL' : indL,
                'IndicesR' : indR,
                'Data' : self._data.array,
                'Mesh' : self.mesh,
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
 
    #-------------------------------------------------
                
    def __lazy_expr_eval_context__(self) : 

        class ctx : 
            def __init__ (self, G) : 
                self.G = G
            def __eq__ (self, y) :
                return isinstance(y, self.__class__) and self.G._is_compatible_for_ops(y.G)
            def __call__ (self, x) : 
                if not isinstance(x, Descriptors.base) : return x
                tmp = self.G.copy()
                x(tmp)
                return tmp
        
        return ctx(self)

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
                if not isinstance(x, Descriptors.base) : return x
                tmp = self.copy()
                x(tmp)
                return tmp
            #e_t2 = self.__lazy_expr_eval_context__()
            self.copyFrom ( lazy_expressions.eval_lazy_expr(e_t, A2) )
        elif isinstance(A, lazy_expressions.lazy_expr_terminal) : #e.g. g<<= SemiCircular (...) 
            self <<= lazy_expressions.lazy_expr(A)
        elif Descriptors.is_scalar(A) : #in the case it is a scalar .... 
            self <<= lazy_expressions.lazy_expr(A)
        elif isinstance(A, GF_Initializers.base) : # backwards compatibility, deprecated
            A(self)
        else :
            raise RuntimeError, " GF Block : <<= operator : RHS not understood"
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
            c = self.copy_with_new_stat(GF_Statistic.Boson if self.mesh.Statistic == y.mesh.Statistic else GF_Statistic.Fermion)
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

    def replaceByTail(self,start) : 
        d = self._data.array
        t = self._tail
        for n, om in enumerate(self.mesh) : # not the most efficient ...
            if n >= start : d[:,:,n] = t(om).array

    #-------------------------------------------------------------------

    # remove this
    def copy_and_transpose(self):
        """deprecated :Transposes the GF Bloc, but gives back a new instance. Does not change the actual GF!"""
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
        """Transposes the GF Bloc: return a new transposed view"""
        
        return self.__class__( 
                Indices = list(self.Indices),
                Mesh  = self.mesh,
                Data = self._data.array.transpose( (1,0,2) ), 
                Tail = self._tail.transpose(),
                Name = self.Name+'(t)', 
                Note = self.Note)

    def conjugate(self):
        """Complex conjugate of the GF Bloc. It follow the policy of numpy and
        make a copy only if the Green function is complex valued"""
        
        return self.__class__( 
                Indices = list(self.Indices),
                Mesh  = self.mesh,
                Data = self._data.array.conjugate(), 
                Tail = self._tail.conjugate(self.mesh.TypeGF==GF_Type.Imaginary_Frequency),
                Name = self.Name+'*', 
                Note = self.Note)

    #-----------------------------------------------------

    # Put it out a a free funciton in a module. Nothing to do here...
    # why is this here ???
    def Delta(self) :
        """Computes Delta from self ...."""
        if self.mesh.TypeGF not in [GF_Type.Real_Frequency, GF_Type.Imaginary_Frequency] :
            raise RuntimeError, "Can not compute Delta for this GF"
        G0 = self if self._tail.OrderMin <=-1 else inverse(self)
        tmp = G0.copy()
        tmp <<= GF_Initializers.A_Omega_Plus_B(G0._tail[-1], G0._tail[0])
        tmp -= G0
        return tmp

#-----------------------------------------------------

from pytriqs.base.Archive.HDF_Archive_Schemes import register_class
register_class (TailGF)
register_class (MeshGF)


