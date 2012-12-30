
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
from ArrayViewWithIndexConverter import ArrayViewWithIndexConverter,_IndicesConverter
import lazy_expressions, Descriptors
from pytriqs.base.Utility.myUtils import *
from pytriqs.base.Plot.protocol import clip_array

class _Plot_Wrapper_Partial_Reduce : 
    """ Internal Use"""
    def __init__(self, obj,  **opt) : 
        self.obj, self.opt = obj,opt
    def _plot_(self,Options) : 
        Options.update(self.opt)
        return self.obj._plot_(Options)

class _GFBloc_concept_impl:

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
    
    Indices = property(__get_Indices)
    IndicesL = property(__get_IndicesL)
    IndicesR = property(__get_IndicesR)
            
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
	return self.Name if self.Name else repr(self)
    
    #-------------------------------------------------

    def __repr__(self) : 
	return """%s %s :  Beta = %.3f; IndicesL = %s, IndicesR = %s """%(self.__class__.__name__, self.Name,
          self.Beta, [x for x in self.IndicesL], [x for x in self.IndicesR])

    #-------------------------------------------------

    def __iter__(self) :
	for i in self._IndicesL : 
	    for j in self._IndicesR :
                b =self[i,j]
		b.Name = "%s_%s_%s"%(self.Name if hasattr(self,'Name') else '',i,j)
		yield i,j,b

    #-----------------------------------------------------
 
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
 
        
