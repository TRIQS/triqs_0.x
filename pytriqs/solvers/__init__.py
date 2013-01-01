
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

__add__ =[] 

from types import *
from pytriqs.base.utility.my_utils import *

class SolverBase(object, Pretty_Print):
    """
        - GFstruct : structure of the Green function in the format
                         (  ( name_of_the_bloc , [list of names for the bloc components]), ...)
                        e.g.( ("up" , [1,2,3,4]), ("down", [1,2,3,4]) ) 
    """

    def __init__(self,GFstruct,param):
        """ GFstruct : the Green function structure as list of tuple"""
        from copy import deepcopy
        self.__param_at_init = deepcopy(param)
        self.__dict__.update(param)
        self.GFStruct = GFstruct[:]
        
        # Name grabbed by the solver
        if type(self.GFStruct[0][1][0])==type(()):
          self.GF_Names   = [ (cle, ["%s-%s"%s for s in val])  for cle,val in self.GFStruct]
        elif [s for s in [1,'1'] if type(self.GFStruct[0][1][0])==type(s) ]:
          self.GF_Names   = [ (cle, ["%s"%s for s in val])  for cle,val in self.GFStruct]
        else: raise TypeError,"Solver indices must be a tuple, int, or string"
        self.Converged = False
        self.Name =''

    def need_special_last_iter(self) : return False

