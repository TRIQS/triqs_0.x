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
import GF_Initializers
import lazy_expressions,Descriptors
from types import IntType,SliceType,StringType
from tools import PlotWrapperPartialReduce, lazy_ctx, IndicesConverter

def py_deserialize( cls, s) : 
    return cls(boost_serialization_string = s)

# Function that transcribe the indices to C++
cdef indices_2_t make_c_indices(obj) : 
    cdef vector[vector[std_string]] res
    cdef vector[std_string] vl,vr
    for i in obj.indicesL : vl.push_back(i)
    for i in obj.indicesR : vr.push_back(i)
    res.push_back(vl); res.push_back(vr)
    return indices_2_t(res)

cdef class _ImplGfLocal :
    cdef object _myIndicesGFBlocL, _myIndicesGFBlocR, _Name, dtype, _IndicesR, _IndicesL,__indices_converter
    def __init__(self, d) : 
        
        # exclusive : size = (n1,n2) or IndicesL/R
        self._IndicesL = list ( d.pop('IndicesL',()) or d.pop('Indices',()) )
        self._IndicesR = list ( d.pop('IndicesR',()) or self._IndicesL  )
        self._Name = d.pop('Name','g')

        # Now check the indices
        ty = set([type(x) for x in self._IndicesL]+[type(x) for x in self._IndicesR])
        assert len(ty)==1, " All indices must have the same type"

        #if ty.pop() != StringType : 
        #    assert self._IndicesL == range (len(self._IndicesL)), " Indices are string, or a range of integer"

        # If the indices are not string, make them string anyway
        self._IndicesL = [ str(x) for x in self._IndicesL ]     
        self._IndicesR = [ str(x) for x in self._IndicesR ]     

        self.__indices_converter = [ IndicesConverter(self._IndicesL), IndicesConverter(self._IndicesR)]

    #-------------  Indices management ---------------------

    property Name : 
        """Name of the Green function (for plots, etc...) """
        def __get__(self) : return self._Name
        def __set__(self,val) : self._Name = str(val)
    
    property indices : 
        """Indices ..."""
        def __get__(self) : 
            if self._IndicesL != self._IndicesR : raise RuntimeError, "Indices R and L are not the same. I can not give you the Indices"
            for ind in self._IndicesL: 
                yield ind

    property indicesL : 
        """Indices ..."""
        def __get__(self) : 
            for ind in self._IndicesL: 
                yield ind
    
    property indicesR : 
        """Indices ..."""
        def __get__(self) : 
            for ind in self._IndicesR: 
                yield ind

    #---------------------   [  ] operator        ------------------------------------------
    
    def __getitem__(self,key):
        """Key is a tuple of index (n1,n2) as defined at construction"""
        if len(key) !=2 : raise KeyError, "[ ] must be given two arguments"
        sl1, sl2 = key
        if type(sl1) == StringType and type(sl2) == StringType :
            # Convert the indices to integer
            sl1, sl2 =  [ self.__indices_converter[i].convertToNumpyIndex(k) for i,k in enumerate(key) ]
        return self.__class__(IndicesL = self._IndicesL[sl1],
                              IndicesR = self._IndicesR[sl2],
                              Name = self.Name,
                              Mesh = self.mesh,
                              Data = self.data[sl1,sl2,:],
                              Tail = self.tail._make_slice(sl1,sl2))

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
 
    #---------------- Repr, str ---------------------------------

    def __str__ (self) : 
        return self.Name if self.Name else repr(self)
    
    def __repr__(self) : 
        return """%s %s : IndicesL = %s, IndicesR = %s"""%(self.__class__.__name__, self.Name,
          [x for x in self._IndicesL], [x for x in self._IndicesR])

    #--------------   PLOT   ---------------------------------------

    property real : 
        """Use self.real in a plot to plot only the real part"""
        def __get__ (self): return PlotWrapperPartialReduce(self,RI='R')

    property imag : 
        """Use self.imag in a plot to plot only the imag part"""
        def __get__ (self): return PlotWrapperPartialReduce(self,RI='I')
  
    def _plot_base (self, OptionsDict, xlabel, ylabel, use_ris, X):
        """ Plot protocol. OptionsDict can contain : 
             * :param RIS: 'R', 'I', 'S', 'RI' [ default] 
             * :param x_window: (xmin,xmax) or None [default]
             * :param Name: a string [default ='']. If not '', it remplaces the name of the function just for this plot.
        """
        import impl_plot
        impl_plot.plot_base(self,OptionsDict,  xlabel, ylabel, use_ris, X)

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
        elif isinstance(A, GF_Initializers.Base) : # backwards compatibility, deprecated
            A(self)
        else :
            raise RuntimeError, " GF Block : <<= operator : RHS not understood"
        return self

    #---------------------------------------------------

    def transpose(self):
        """Transposes the GF Bloc: return a new transposed view"""
       
       ### WARNING : this depends on the C++ layering ....
        return self.__class__( 
                Indices = list(self.Indices),
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
                Indices = list(self.Indices),
                Mesh  = self.mesh,
                Data = self.data.conjugate(), 
                Tail = self.tail.conjugate(), #self.mesh.TypeGF==GF_Type.Imaginary_Frequency),
                Name = self.Name+'*', 
                Note = self.Note)
    
    #------------------  Density -----------------------------------
    
    def total_density(self) :
        """Trace density"""
        return numpy.trace(self.density())

