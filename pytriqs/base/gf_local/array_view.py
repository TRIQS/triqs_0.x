
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

import numpy
from types import *

class _IndicesConverter : 
    def __init__(self, indices) : 
        if indices == None : # none is a trivial converter
            self.indices=None
            return

        #print "Indices are ", indices
        try : 
            self.indices = list(indices)[:] # make a copy
        except :
            raise RuntimeError, "Indices must be a list or transformable into a list %s"(indices,)
        assert len(set(repr(x) for x in self.indices)) == len(self.indices), "Error : indices are not unique !!!"
        assert self.indices != [], "Error : indices are empty!!!"
        try :
            # a continuous list of ordered integers
            self.__indices_are_sliceable = (self.indices== range(min(self.indices),max(self.indices)+1)) 
            self.__indexmin = min(self.indices)
            self.__indexmax = max(self.indices)
            #self.__indexlen = self.__indexmax - self.__indexmin +1
        except:
            self.__indices_are_sliceable =False
        self.Length = len(self.indices)
    
    def convertToNumpyIndex(self,a,noslice=False):
        """
        Transcription of one Python index/slice into the index/slice of the numpy
        """
        if self.indices==None: 
            if type(a)==SliceType : return a
            return a if noslice else slice(a,a+1,1)

        if type(a)==SliceType :
            if not self.__indices_are_sliceable : raise IndexError, "Indices %s can't be sliced"%self.indices
            # put the slice to start at 0, and clip it
            sta,sto,ste = slice(a.start, a.stop , a.step).indices(self.__indexmax+1)
            return slice( max(0,sta-self.__indexmin), sto - self.__indexmin, ste)
                    
        try :
            s1 = self.indices.index(a) 
        except  :
            raise IndexError, "Index %s out of range %s"%(a,self.indices)
        return s1 if noslice else slice(s1,s1+1,1)

######################################################
   
class ArrayViewWithIndexConverter(object):
    """
    This class is a view of a numpy array, with a converter for the indices.
    Construction :
       ArrayViewWithIndexConverter ( A , indices)
       where indices are an index generator, and A is an array.

    Usage : 
       # g is a GFBloc with indices [1,2,...] or ['A','B', ...]
       
       # Create such an array from g  e.g.
       A = g.array_with_indices()
       A['a','b'] = 2

       # When you slice, you get back a ArrayViewWithIndexConverter: 
       type(A[:,1])

       NB : This can be quite slower that direct calculation of the array, because of the indices conversions.
           In that case, use A.array (readonly property) which returns the underlying numpy array.
           (in the correct order).
    """
    def __init__(self, A, *IndicesList):
        """ 
        Inputs : 
          * A is the array to be viewed
          * IndicesList is a list of names of indices for first, second, ... indices
          """
        self.__converters = [ _IndicesConverter(ind) for ind in IndicesList]
        self.__A = A 
        self.__dim= len(A.shape)
        assert self.__dim>= len( IndicesList), "Error : you provided to many indices for this array !"
        #assert self.__dim<= len( IndicesList), "Error : you provided to few indices for this array !"

        for n,conv in enumerate(self.__converters) : 
            if conv.indices and conv.Length !=  A.shape[n] :
                print conv.indices,A
                raise RuntimeError, "Array size and number of indices mismatch for index %s"%n
        
    def __getitem__(self,key):
        assert len(key)==self.__dim , "Indices must be of the form (n1,n2,...) instead of %s"%(key,)
        sl,ind = [],[]
        for (c,k) in zip(self.__converters,key) :
            s = c.convertToNumpyIndex(k,noslice=True) #either an index or a slice
            sl.append(s)
            if type(s) ==SliceType :
                ind.append(c.indices[s] if c.indices else None)
        #print "ind is ",ind,sl
        return ArrayViewWithIndexConverter( self.__A[tuple(sl)], *ind) if ind else self.__A[tuple(sl)]
        
    def __setitem__(self,key,val):
        assert len(key)==self.__dim , "Indices must be of the form (n1,n2,...) instead of %s"%(key,)
        sl = tuple([ c.convertToNumpyIndex(k) for (c,k) in zip(self.__converters,key)]) # always a slice !
        self.__A[sl]= val
    
    def __get_array(self) :
        return self.__A
    
    array = property(__get_array,None)

    def Indices(self,n) :
        for i in self.__converters[n].indices if self.__converters[n].indices else range(self.__A.shape[n]) : 
            yield i
             
    def __repr__(self) :
        ind = dict([ (n,c.indices if c.indices else range(self.__A.shape[n])) for n,c in enumerate(self.__converters)])
        return "ArrayViewWithIndexConverter with : \n  Indices = %s\n  Array = %s"%(ind,self.__A)
    
    
