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

# Function that transcribe the indices to C++
cdef indices_2_t make_c_indices(indicesL, indicesR) : 
    cdef vector[vector[std_string]] res
    cdef vector[std_string] vl,vr
    for i in indicesL : vl.push_back(i)
    for i in indicesR : vr.push_back(i)
    res.push_back(vl); res.push_back(vr)
    return indices_2_t(res)

cdef class GfGeneric_cython :
    
    cdef object _name, dtype, _mesh, _data, _singularity, _symmetry, _indices
    cdef readonly  _derived

    def __init__(self, mesh, data, singularity, symmetry, indices, name, derived ) : 
        self._mesh, self._data, self._singularity, self._symmetry, self._indices= mesh, data, singularity, symmetry, indices
        self._name = name
        self._derived = derived

    def __reduce_to_dict__(self):
        return { 'Mesh' : self._mesh, 'Data' : self._data, 
                'Tail' : self._singularity, 'Symmetry' : self._symmetry,
                'Indices' : self._indices, 'Name' : self._name } 

    property mesh : 
        """Mesh"""
        def __get__(self): return self._mesh
    
    property tail : 
        def __get__(self): return self._singularity
        def __set__(self,TailGf t): 
            assert (self.N1, self.N2, self._singularity.size) == (t.N1, t.N2, t.size)
            self._singularity.copy_from (t)

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

    property name : 
        """Name of the Green function (for plots, etc...) """
        def __get__(self) : return self._name
        def __set__(self,val) : self._name = str(val)

    property Name : 
        """Name of the Green function (for plots, etc...) """
        # DEPRECATED
        def __get__(self) : return self._name

