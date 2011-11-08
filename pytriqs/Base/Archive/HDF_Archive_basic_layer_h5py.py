
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

# h5py
import numpy,string,h5py
class HDF_Archive_group_basic_layer : 
    _class_version = 1

    def __init__(self, parent, subpath ): 
        """  """
        self.options = parent.options
        self._group = parent._group[subpath] if subpath else parent._group
        assert type(self._group) in [h5py.highlevel.Group,h5py.highlevel.File], "Internal error"
        self.ignored_keys = [] 

    def _init_root(self, LocalFileName, Open_Flag) : 
        try : 
            fich = h5py.File(LocalFileName, Open_Flag)
        except : 
            print "Can not open the HDF file %s"%LocalFileName
            raise
        # checking the version
        if Open_Flag not in ['r','r+','a'] : 
            self._version = self._class_version
        else : 
            try : 
                self._version = int(fich.attrs['HDF_Archive_Version']) 
            except : 
                self._version = 1
            if self._version > self._class_version : 
                raise IOError, "File %s is too recent for this version of HDF_Archive module"%Filename
        self._group = fich

    def is_group(self,p) :
        """Is p a subgroup ?"""
        assert len(p)>0 and p[0]!='/'
        return p in self._group and type(self._group[p]) == h5py.highlevel.Group
    
    def is_data(self,p) :
        """Is p a leaf ?"""
        assert len(p)>0 and p[0]!='/' 
        return p in self._group and type(self._group[p]) == h5py.highlevel.Dataset

    def write_attr (self, key, val) : 
        self._group.attrs[key] =  val 

    def read_attr(self,AttributeName) : 
        return self._group.attrs[AttributeName] 

    def _read (self, key) : 
        A = self._group[key] 
        val =  numpy.array(A) if A.shape!=() else A.value 
        if self.options["UseAlpsNotationForComplex"] and '__complex__' in self._group[key].attrs :
            assert type(val) == numpy.ndarray, 'complex tag is set, but I have not an array'
            assert not numpy.iscomplexobj(val), 'complex tag is set, but I have a complex !'
            if len(val.shape) == 1 : 
                val = val[0] + 1j * val[1]
            else : 
                val = val[...,0] + 1j*val[...,1]

        def _numpy_scalar_to_python_scalar (v) :
            n= numpy
            if type(v) in [n.int, n.int8, n.int16, n.int32, n.int64] : return int(v)
            if type(v) in [n.uint8, n.uint16, n.uint32, n.uint64] : return int(v)
            if type(v) in [n.float,  n.float32, n.float64] : return float(v)
            if type(v) in [n.complex, n.complex64, n.complex128] : return complex(v)
            if type(v) in [n.str_, n.string_] : return str(v)
            return v 
  
        return _numpy_scalar_to_python_scalar ( val) 
      
    def _write_array(self, key, A) :
        c =  self.options["UseAlpsNotationForComplex"] and numpy.iscomplexobj(A)  
        if c: 
            val  = numpy.zeros( A.shape + (2,) ) 
            val[...,0], val[...,1] = A.real, A.imag
        else : 
            val = A
        self._group[key] = numpy.array(val,copy=1,order='C')  
        if c : self._group[key].attrs["__complex__"] = 1

    def _write_scalar(self, key, A) :
        c =  self.options["UseAlpsNotationForComplex"] and type(A) ==type (1j) 
        val = numpy.array([A.real, A.imag]) if c else A
        self._group[key] =val 
        if c : self._group[key].attrs["__complex__"]= 1
  
    def _flush(self) : 
        self._group.file.flush()
  
    def create_group (self,key):
        self._group.create_group(key)

    def _keys(self) :
        def res() : 
            for name in  self._group.iterkeys():
                yield name
        return res()

    def _clean_key(self,key, report_error=False) :
        if report_error and key not in self._group : 
             raise KeyError, "Key %s is not in archive !!"%key 
        if key in self._group : del self._group[key]
        else: raise KeyError, "Key %s is not in archive !!"%key 
   
