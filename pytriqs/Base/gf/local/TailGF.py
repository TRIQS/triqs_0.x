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
from pytriqs_GF import TailGF
import numpy
from math import pi
from pytriqs.Base.GF_Local.ArrayViewWithIndexConverter import ArrayViewWithIndexConverter

#-----------------------------------------------------
#  Code Injection
#-----------------------------------------------------

from pytriqs.Base.Utility.Injector import make_injector        # inject new code in the SAME class  

class __inject (make_injector(TailGF), TailGF):
    """ 
    """
    
    def __init__(self, **d):
        """
        """
        self._indL, self._indR = d['IndicesL'], d['IndicesR']
        N1,N2 = len(self._indL), len(self._indR)
        self.__data_array = d['array'] if 'array' in d else numpy.zero((N1,N2,d['size'])) 
        self._omin = d['OrderMin']
        self._init_before_injection__(self.__data_array, self._omin) 

    def copy(self) : 
        return self.__class__(IndicesL = self._indL, IndicesR = self._indR, array = self.__data_array.copy(), OrderMin = self._omin)
 
    def __repr__ (self) :
        return string.join([ "%s"%self[r-self._omin] + ("/" if r<0 else "") + "Om^%s"%(abs(r)) for r in range (self._omin, self._omax+1) ] , "+") 

    def __getitem__(self,i) :
        """Returns the i-th coefficient of the expansion, or order Om^i"""
        if not self.has_coef(i) : raise IndexError, "Index %s is out of range"%i
        return ArrayViewWithIndexConverter(self._data.array[ i - self._omin, :,:], self.IndicesL, self.IndicesR)

    def __setitem__(self,i, val) :
        """Sets the i-th coefficient of the expansion, or order Om^i"""
       if not self.has_coef(i) : raise IndexError, "Index %s is out of range"%i
       self._data.array[ i - self._omin, :,:] = val

    @property
    def OrderMin(self) : return self._omin

    @property
    def OrderMax(self) : return self._omin+self.size() -1

    @property
    def size(self) : return self.__data_array.shape[2]

    #---------------------------------------------------------------------------------
 
    @property
    def _data(self) : 
        return ArrayViewWithIndexConverter(self.__data_array, self._IndicesL, self._IndicesR, None)

    @_data.setter
    def _data(self, value) : self.__data_array[:,:,:] = value 
 
    #-----------------------------------------------------

    def __reduce_to_dict__(self):
        indL = repr(tuple(self._indL))
        indR = repr(tuple(self._indR))
        assert(eval(indL)==tuple(self._indL))
        assert(eval(indR)==tuple(self._indR))
        return {'IndicesL' : indL,
                'IndicesR' : indR,
                'array' : self._data.array,
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
    def __call__(self, n) : return self.__data_array[:,:,n - self.OrderMin]

    def zero(self) : 
        """Sets the expansion to 0"""
        self.__data_array[:,:,:] =0

    #.def ("copyFrom",&tail_view::copyFrom)
    #.def ("invert",&tail_view::invert,"Replace itself by the expansion of the inverse of the function.")
 
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


