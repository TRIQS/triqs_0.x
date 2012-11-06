
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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


import copy,numpy
import string
from types import *
from pytriqs.Base.GF_Local.GF import GF
from pytriqs.Base.Archive.HDF_Archive import *
import pytriqs.Base.Utility.MPI as MPI


class Symmetry:
    """This class provides the routines for applying symmetry operations for the k sums.
       It contains the permutations of the atoms in the unti cell, and the corresponding
       rotational matrices for each symmetry operation."""

    def __init__(self, HDFfile, subgroup = None):
        """Initialises the class.
           Reads the permutations and rotation matrizes from the file, and constructs the mapping for
           the given orbitals. For each orbit a matrix is read!!!
           SO: Flag for SO coupled calculations.
           SP: Spin polarisation yes/no
           """

        assert type(HDFfile)==StringType,"HDFfile must be a filename"; self.HDFfile = HDFfile
        thingstoread = ['Ns','Natoms','perm','orbits','SO','SP','timeinv','mat','mat_tinv']
        for it in thingstoread: exec "self.%s = 0"%it

        if (MPI.IS_MASTER_NODE()):
            #Read the stuff on master:
            ar = HDF_Archive(HDFfile,'a')
            if (subgroup is None):
                ar2 = ar
            else:
                ar2 = ar[subgroup]

            for it in thingstoread: exec "self.%s = ar2['%s']"%(it,it)
            del ar2
            del ar

        #broadcasting
        for it in thingstoread: exec "self.%s = MPI.bcast(self.%s)"%(it,it)
        
        # now define the mapping of orbitals:
        # self.map[iorb]=jorb gives the permutation of the orbitals as given in the list, when the 
        # permutation of the atoms is done:
        self.N_orbits = len(self.orbits)

        self.map = [ [0 for iorb in range(self.N_orbits)] for iNs in range(self.Ns) ]
        for iNs in range(self.Ns):
            for iorb in range(self.N_orbits):
             
                srch = copy.deepcopy(self.orbits[iorb])
                srch[0] = self.perm[iNs][self.orbits[iorb][0]-1]
                self.map[iNs][iorb] = self.orbits.index(srch)
                    
       

    def symmetrise(self,obj):
        
        assert isinstance(obj,list),"obj has to be a list of objects!"
        assert len(obj)==self.N_orbits,"obj has to be a list of the same length as defined in the init"

        if (isinstance(obj[0],GF)):
            symm_obj = [ obj[i].copy() for i in range(len(obj)) ]        # here the result is stored, it is a GF!
            for iorb in range(self.N_orbits): symm_obj[iorb].zero()      # set to zero
        else:
            # if not a GF, we assume it is a matrix (density matrix), has to be complex since self.mat is complex!
            #symm_obj = [ numpy.zeros([self.orbits[iorb][3],self.orbits[iorb][3]],numpy.complex_) for iorb in range(self.N_orbits) ]
            symm_obj = [ copy.deepcopy(obj[i]) for i in range(len(obj)) ]
         
            for iorb in range(self.N_orbits):
                if (type(symm_obj[iorb])==DictType):
                    for ii in symm_obj[iorb]: symm_obj[iorb][ii] *= 0.0
                else:
                    symm_obj[iorb] *= 0.0
        
                
        for iNs in range(self.Ns):

            for iorb in range(self.N_orbits):

                l = self.orbits[iorb][2]         # s, p, d, or f
                dim = self.orbits[iorb][3]
                jorb = self.map[iNs][iorb]

             
                if (isinstance(obj[0],GF)):
                    
                    #if l==0:
                    #    symm_obj[jorb] += obj[iorb]
                    #else:
                    
                    tmp = obj[iorb].copy()
                    if (self.timeinv[iNs]): tmp <<= tmp.transpose()
                    for sig,gf in tmp: tmp[sig].from_L_G_R(self.mat[iNs][iorb],tmp[sig],self.mat[iNs][iorb].conjugate().transpose())
                    tmp *= 1.0/self.Ns
                    symm_obj[jorb] += tmp

                else:

                    if (type(obj[iorb])==DictType):

                        for ii in obj[iorb]:
                            #if (l==0):
                            #    symm_obj[jorb][ii] += obj[iorb][ii]/self.Ns
                            #else:
                            if (self.timeinv[iNs]==0):
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[iNs][iorb],obj[iorb][ii]),
                                                                self.mat[iNs][iorb].conjugate().transpose()) / self.Ns
                            else:
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[iNs][iorb],obj[iorb][ii].conjugate()),
                                                                self.mat[iNs][iorb].conjugate().transpose()) / self.Ns

                            

                    else:
                        #if (l==0):
                        #    symm_obj[jorb] += obj[iorb]/self.Ns
                        #else:
                        if (self.timeinv[iNs]==0):
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[iNs][iorb],obj[iorb]),self.mat[iNs][iorb].conjugate().transpose()) / self.Ns
                        else:
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[iNs][iorb],obj[iorb].conjugate()),
                                                        self.mat[iNs][iorb].conjugate().transpose()) / self.Ns
                        
     
# This does not what it is supposed to do, check how this should work:       
#        if ((self.SO==0) and (self.SP==0)):
#            # add time inv:
            #MPI.report("Add time inversion")
#            for iorb in range(self.N_orbits):
#                if (isinstance(symm_obj[0],GF)):
#                    tmp = symm_obj[iorb].copy()
#                    tmp <<= tmp.transpose()
#                    for sig,gf in tmp: tmp[sig].from_L_G_R(self.mat_tinv[iorb],tmp[sig],self.mat_tinv[iorb].transpose().conjugate())
#                    symm_obj[iorb] += tmp
#                    symm_obj[iorb] /= 2.0
#                    
#                else:
#                    if (type(symm_obj[iorb])==DictType):
#                        for ii in symm_obj[iorb]:
#                            symm_obj[iorb][ii] += numpy.dot(numpy.dot(self.mat_tinv[iorb],symm_obj[iorb][ii].conjugate()),
#                                                            self.mat_tinv[iorb].transpose().conjugate())
#                            symm_obj[iorb][ii] /= 2.0
#                    else:
#                        symm_obj[iorb] += numpy.dot(numpy.dot(self.mat_tinv[iorb],symm_obj[iorb].conjugate()),
#                                                    self.mat_tinv[iorb].transpose().conjugate())
#                        symm_obj[iorb] /= 2.0
                                
    
        return symm_obj
                    
  
                    
 
