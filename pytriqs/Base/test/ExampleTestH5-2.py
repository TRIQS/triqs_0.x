
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

from pytriqs.Base.Archive import *
from numpy import *


class test_p2: 
    """ This class creates a group"""
    def __init__(self,l=None) :
        if l : 
          self.a = l*array([[1,2,3],[4,5,6]])
          self.b = "un test"

    def __str__(self): 
        return """
        a = %s
        b = %s"""%(self.a,self.b)

    def __write_hdf5__(self,archive, path) : 
        archive.create_group(path)
        archive.write("%s/a"%path,self.a)
        archive.write("%s/b"%path,self.b)

    def __read_hdf5__(self,archive, path) : 
        if not archive.is_group(path) : raise IOError, "Pb"
        self.a = archive.read ("%s/a"%path)
        self.b = archive.read ("%s/b"%path)

class test_p2B: 
    """This class does not create a group"""
    def __init__(self,l=None) :
        if l : 
          self.b = l

    def __str__(self): 
        return """
        b = %s"""%(self.a,self.b)

    def __write_hdf5__(self,archive, path) : 
        archive.write("%s"%path,self.b)

    def __read_hdf5__(self,archive, path) : 
        self.b = archive.read ("%s"%path)

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (test_p2)
register_class (test_p2B)

h = HDF_Archive('ExampleTestH5-2.output.h5','w')

t1 = test_p2(1)
t2 = test_p2(2)
h['t'] = t1
h['t'] = t2
h['B'] = test_p2B(123)

h = HDF_Archive('ExampleTestH5-2.output.h5','r')
print h['t']
print h['B']


