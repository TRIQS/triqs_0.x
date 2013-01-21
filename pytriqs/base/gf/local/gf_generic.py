
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
import lazy_expressions, descriptors
import pytriqs.base.utility.my_utils
from pytriqs.base.plot.protocol import clip_array
from types import IntType, SliceType, StringType
from tools import LazyCTX, IndicesConverter, get_indices_in_dict, py_deserialize
from impl_plot import PlotWrapperPartialReduce

class GfGeneric:

    def copy(self):
        return self._derived(indices = self.indices, mesh = self.mesh, data = self.data.copy(), tail = self.tail.copy(), name = self.name)

    def copy_from(self, X):
        assert self._derived is X._derived
        assert self.mesh == X.mesh
        self.data = X.data
        self.tail = X.tail
        #assert list(self._indices)== list(X._indices)
        #self._symmetry = X._symmetry
        self._name = X.name

    #---------------------   [  ] operator        ------------------------------------------

    def __getitem__(self, key):
        """Key is a tuple of index (n1, n2) as defined at construction"""
        if len(key) !=2: raise KeyError, "[ ] must be given two arguments"
        sl1, sl2 = key
        if type(sl1) == StringType and type(sl2) == StringType:
            # Convert the indices to integer
            indices_converter = [ IndicesConverter(self.indicesL), IndicesConverter(self.indicesR)]
            sl1, sl2 =  [ indices_converter[i].convertToNumpyIndex(k) for i, k in enumerate(key) ]
        if type (sl1) != slice: sl1 = slice (sl1, sl1+1)
        if type (sl2) != slice: sl2 = slice (sl2, sl2+1)
        return self.__class__(indicesL = list(self.indicesL)[sl1],
                              indicesR = list(self.indicesR)[sl2],
                              name = self.name,
                              mesh = self.mesh,
                              data = self.data[sl1, sl2,:],
                              tail = self.tail._make_slice(sl1, sl2))

    def __setitem__(self, key, val):
        g = self.__getitem__(key)
        g <<= val

    #------------- Iteration ------------------------------------

    def __iter__(self):
        for i in self.indicesL:
            for j in self.indicesR:
                b =self[i, j]
                b.name = "%s_%s_%s"%(self.name if hasattr(self, 'name') else '', i, j)
                yield i, j, b

    #---------------- Repr, str ---------------------------------

    def __str__(self):
        return self.name if self.name else repr(self)

    def __repr__(self):
        return """%s %s: indicesL = %s, indicesR = %s"""%(self.__class__.__name__, self.name,
          [x for x in self.indicesL], [x for x in self.indicesR])

    #--------------   PLOT   ---------------------------------------

    @property
    def real(self):
        """Use self.real in a plot to plot only the real part"""
        return PlotWrapperPartialReduce(self, RI='R')

    @property
    def imag(self):
        """Use self.imag in a plot to plot only the imag part"""
        return PlotWrapperPartialReduce(self, RI='I')

    #------------------

    def x_data_view(self, x_window = None, flatten_y = False):
        """
        :param x_window: the window of x variable (omega/omega_n/t/tau) for which data is requested
                          if None, take the full window
        :param flatten_y: If the Green function is of size (1, 1) flatten the array as a 1d array
        :rtype: a tuple (X, data) where
                 * X is a 1d numpy of the x variable inside the window requested
                 * data is a 3d numpy array of dim (:,:, len(X)), the corresponding slice of data
                   If flatten_y is True and dim is (1, 1, *), returns a 1d numpy
        """
        X = [x.imag for x in self.mesh] if type(self.mesh) == MeshImFreq else [x for x in self.mesh]
        X, data = numpy.array(X), self.data
        if x_window:
          sl = clip_array (X, *x_window) if x_window else slice(len(X)) # the slice due to clip option x_window
          X, data = X[sl],  data[:,:, sl]
        if flatten_y and data.shape[:2]==(1, 1): data = data[0, 0,:]
        return X, data

    #--------  LAZY expression system -----------------------------------------

    def __lazy_expr_eval_context__(self):
        return LazyCTX(self)

    def __eq__(self, other):
        raise RuntimeError, " Operator not defined "

    def __ilshift__(self, A):
        """ A can be two things:
          * G <<= any_gf_init will init the GFBloc with the initializer
          * G <<= g2 where g2 is a GFBloc will copy g2 into self
        """
        if isinstance(A, self.__class__):
            if self is not A: self.copy_from(A) # otherwise it is useless AND does not work !!
        elif isinstance(A, lazy_expressions.LazyExpr): # A is a lazy_expression made of GF, scalars, descriptors
            A2= descriptors.convert_scalar_to_const(A)
            def e_t (x):
                if not isinstance(x, descriptors.Base): return x
                tmp = self.copy()
                x(tmp)
                return tmp
            self.copy_from (lazy_expressions.eval_expr_with_context(e_t, A2) )
        elif isinstance(A, lazy_expressions.LazyExprTerminal): #e.g. g<<= SemiCircular (...)
            self <<= lazy_expressions.LazyExpr(A)
        elif descriptors.is_scalar(A): #in the case it is a scalar ....
            self <<= lazy_expressions.LazyExpr(A)
        else:
            raise RuntimeError, " <<= operator: RHS  not understood"
        return self

    #--------------------  Arithmetic operations  ---------------------------------

    def __iadd__(self, arg):
        d, t = self.data, self.tail
        if type(self) == type(arg):
            d[:,:,:] += arg.data
            t += arg.tail
        elif isinstance(arg, numpy.ndarray): # an array considered as a constant function
            for om in range (d.shape[-1]): d[:,:, om ] += arg
            t[0][:,:] += arg
        elif descriptors.is_scalar(arg): # just a scalar
            arg = arg*numpy.identity(self.N1)
            for om in range (d.shape[-1]): d[:,:, om ] += arg
            t[0][:,:] += arg
        else:
            raise RuntimeError, " argument type not recognized in += for %s"%arg
        return self

    def __add__(self, y):
        if descriptors.is_lazy(y): return lazy_expressions.make_lazy(self) + y
        c = self.copy()
        c += y
        return c

    def __radd__(self, y): return self.__add__(y)

    def __isub__(self, arg):
        d, t = self.data, self.tail
        if type(self) == type(arg):
            d[:,:,:] -= arg.data
            t -= arg.tail
        elif isinstance(arg, numpy.ndarray): # an array considered as a constant function
            for om in range (d.shape[-1]): d[:,:, om ] -= arg
            t[0][:,:] -= arg
        elif descriptors.is_scalar(arg): # just a scalar
            arg = arg*numpy.identity(self.N1)
            for om in range (d.shape[-1]): d[:,:, om ] -= arg
            t[0][:,:] -= arg
        else:
            raise RuntimeError, " argument type not recognized in -= for %s"%arg
        return self

    def __sub__(self, y):
        if descriptors.is_lazy(y): return lazy_expressions.make_lazy(self) - y
        c = self.copy()
        c -= y
        return c

    def __rsub__(self, y):
        c = (-1)*self.copy()
        c += y
        return c

    def __imul__(self, arg):
        """ If arg is a scalar, simple scalar multiplication
            If arg is a BlockGf(any object with data and tail as in GF), they it is a matrix multiplication, slice by slice
        """
        if type(self) == type(arg):
            d, d2 = self.data, arg.data
            assert d.shape == d2.shape, " Green function block multiplication with arrays of different size !"
            for om in range (d.shape[-1]):
                d[:,:, om ] = numpy.dot(d[:,:, om], d2[:,:, om])
            self.tail = arg.tail * self.tail
        elif type(arg).__name__ in ['float', 'int', 'complex']:
            self.data *= arg
            self.tail *= arg
        else:
            raise RuntimeError, " argument type not recognized in *= for %s"%arg
        return self

    def __mul__(self, arg):
        if descriptors.is_lazy(arg): return lazy_expressions.make_lazy(self) * arg
        res = self.copy()
        res *= arg
        return res

    def __rmul__(self, x):
        assert type(x).__name__ in ['float', 'int', 'complex'], "lhs must be a scalar but I found %s"%x
        return self.__mul__(x)

    def imatmul_L(self, L):
        dot = numpy.dot
        assert type(L).__name__ in  ['ndarray', 'matrix']
        A = self.data
        for i in range(A.shape[-1]):
            A[:,:, i] = dot(L, A[:,:, i])
        return self

    def imatmul_R(self, R):
        dot = numpy.dot
        assert type(R).__name__ in  ['ndarray', 'matrix']
        A = self.data
        for i in range(A.shape[-1]):
            A[:,:, i] = dot(A[:,:, i], R)
        return self

    def imatmul_LR(self, L, R):
        dot = numpy.dot
        assert type(R).__name__ in  ['ndarray', 'matrix']
        assert type(L).__name__ in  ['ndarray', 'matrix']
        A = self.data
        for i in range(A.shape[-1]):
            A[:,:, i] = dot(L, dot(A[:,:, i], R))
        return self

    # RENAME THIS !
    def from_L_G_R(self, L, G, R):
        """ For all argument, replace the matrix by L *matrix * R"""
        d, dg = self.data, G.data
        for om in range (d.shape[-1]):
            d[:,:,om] = numpy.dot(numpy.dot(L, dg[:,:,om]), R)
        for order in range(G.tail.data.shape[-1]):
          self.tail.data[:,:,order] = numpy.dot(numpy.dot(L, G.tail.data[:,:,order]), R)

    def __idiv__(self, arg):
        """ If arg is a scalar, simple scalar multiplication
        """
        if descriptors.is_lazy(arg): return lazy_expressions.make_lazy(self) / arg
        n = type(arg).__name__
        if n in ['float', 'int', 'complex']:
            self.data /= arg
            self.tail /= arg
        else:
            raise RuntimeError, " argument type not recognized in imul for %s"%arg
        return self

    def __div__(self, arg):
        assert type(arg).__name__ in  ['float', 'int', 'complex'], "Error in /"
        res = self.copy()
        res /= arg
        return res

    #---------------------------------------------------

    def invert(self):
        """Invert the matrix for all arguments"""
        d = self.data
        for om in range (d.shape[-1]):
            d[:,:, om ] = numpy.linalg.inv(d[:,:, om])
        self.tail.invert()

    #---------------------------------------------------

    def transpose(self):
        """Transposes the GF Bloc: return a new transposed view"""
        ### WARNING: this depends on the C++ layering ....
        return self.__class__(
                indices = list(self.indices),
                mesh = self.mesh,
                data = self.data.transpose( (1, 0, 2) ),
                tail = self.tail.transpose(),
                name = self.name+'(t)')

    #---------------------------------------------------

    def conjugate(self):
        """Complex conjugate of the GF Bloc. It follow the policy of numpy and
        make a copy only if the Green function is complex valued"""

        return self.__class__(
                indices = list(self.indices),
                mesh = self.mesh,
                data = self.data.conjugate(),
                tail = self.tail.conjugate(),
                name = self.name+'*')

    #------------------  Density -----------------------------------

    def total_density(self):
        """Trace density"""
        return numpy.trace(self.density())

