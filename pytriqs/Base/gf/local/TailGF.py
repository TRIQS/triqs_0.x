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

__all__ = ['TailGF']
from gf import TailGF_c
import numpy
import string
from pytriqs.Base.GF_Local.ArrayViewWithIndexConverter import ArrayViewWithIndexConverter

class TailGF (TailGF_c):
    """ 
    """
    
    def __init__(self, **d):
        """
        """
        self._indL, self._indR = d['IndicesL'], d['IndicesR']
        N1,N2 = len(self._indL), len(self._indR)
        data_array = d['array'] if 'array' in d else numpy.zeros((N1,N2,d['size'])  ,numpy.complex) #,order='F')
        self._omin = d['OrderMin']
        TailGF_c.__init__(self,data_array,self._omin)

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

    @property
    def OrderMin(self) : return self._omin

    @property
    def OrderMax(self) : return self._omin+self.size-1

    @property
    def size(self) : return self._data_raw.shape[2]

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
 
    #---- arithmetic operations ----
 
    ########################################## TEMPORARY
    def __iadd__(self,arg):
        self._data_raw += arg._data_raw
        return self
    def __isub__(self,arg):
        self._data_raw -= arg._data_raw
        return self
    def __imul__(self,arg):
        self._data_raw *= arg
        return self
    def __idiv__(self,arg):
        self._data_raw /= arg
        return self
    ########################################## TEMPORARY

    def __add__(self,y):
        c= self.copy()
        c+=y
        return c

    def __sub__(self,y):
        c= self.copy()
        c-=y
        return c

    def __mul__(self,y):
        c= self.copy()
        c*=y
        return c

    def __rmul__(self,x): return self.__mul__(x)
    def __radd__(self,x): return self.__add__(x)
    def __rsub__(self,x): return self.__sub__(x)
   
    def __div__(self,y):
        c= self.copy()
        c/=y
        return c

    #---- other operations ----
    def __call__(self, n) :
        if not self.has_coef(n) : raise IndexError, "Index %s is out of range"%n
        return self._data_raw[:,:,n-self.OrderMin]

    def zero(self) : 
        """Sets the expansion to 0"""
        self._data_raw[:,:,:] =0

    def from_L_T_R(self,L,R) : 
      """Multiply Left & right with Matrices :  G <- L* G2 * R"""
      assert(0)

    def transpose (self) : 
        """Transpose the array : new view as in numpy"""
        assert(0)

    def conjugate(self) : 
        """Transpose the array : new view as in numpy"""
        assert(0)
        # Hum, a pb here : shall we return a new object (view) or do it one site
        # ? (BAD !!).

#-----------------------------------------------------
#  Register the class for HDF_Archive
#-----------------------------------------------------

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (TailGF)

