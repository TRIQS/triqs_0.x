
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
from math import cos, sin


def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
	    yield string.atof(x)


class Symmetry:
    """This class provides the routines for applying symmetry operations for the k sums.
       It contains the permutations of the atoms in the unti cell, and the corresponding
       rotational matrices for each symmetry operation."""

    def __init__(self, Filename, orbits, SO = 0):
        """Initialises the class.
           Reads the permutations and rotation matrizes from the file, and constructs the mapping for
           the given orbitals. For each orbit a matrix is read!!!
           SO: Flag for SO coupled calculations.
           """

        assert type(Filename)==StringType,"LDA_file must be a filename"


        self.N_orbits = len(orbits)
        self.orbits = copy.deepcopy(orbits)
        self.SO = SO
                
        
        R = Read_Fortran_File(Filename)

        try:

            self.Ns = int(R.next())           # Number of symmetry operations

            self.Natoms = int(R.next())       # number of atoms involved
            
            self.perm = [ [int(R.next()) for i in xrange(self.Natoms)] for j in xrange(self.Ns) ]    # list of permutations of the atoms

            if (self.SO):
                # phases and time-inversion for SO symmetry operations
                self.phase   = [ R.next() for j in xrange(self.Ns) ]
                self.timeinv = [ int(R.next()) for j in xrange(self.Ns) ]

            self.mat = []
            for iNs in xrange(self.Ns):
                
                self.mat.append( [ numpy.zeros([self.orbits[orb][3], self.orbits[orb][3]],numpy.complex_) for orb in xrange(self.N_orbits) ] )
                for orb in range(self.N_orbits):
                    if (self.orbits[orb][2]==0): # read p, d, and f matrices, s is unity
                        self.mat[iNs][orb][0,0] = 1.0
                    else:
                        for i in xrange(self.orbits[orb][3]):
                            for j in xrange(self.orbits[orb][3]):
                                self.mat[iNs][orb][j,i] = R.next()            # real part
                        for i in xrange(self.orbits[orb][3]):
                            for j in xrange(self.orbits[orb][3]):
                                self.mat[iNs][orb][j,i] += 1j * R.next()      # imaginary part


        except StopIteration : # a more explicit error if the file is corrupted.
	    raise "Symmetry : reading file failed!"
        
        R.close()

        
        # now define the mapping of orbitals:
        # self.map[iorb]=jorb gives the permutation of the orbitals as given in the list, when the 
        # permutation of the atoms is done:

        self.map = [ [0 for iorb in range(self.N_orbits)] for iNs in range(self.Ns) ]
        for iNs in range(self.Ns):
            for iorb in range(self.N_orbits):
             
                srch = copy.deepcopy(self.orbits[iorb])
                srch[0] = self.perm[iNs][self.orbits[iorb][0]-1]
                self.map[iNs][iorb] = self.orbits.index(srch)
                    
       

    def symmetrise_noSO(self,obj):
        """ symmetrisation without SO coupling"""
        
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

                    if l==0:
                        tmp = obj[iorb].copy()
                        tmp *= 1.0/self.Ns
                        symm_obj[jorb] += tmp
                    else:
                    
                        tmp = obj[iorb].copy()
                        for sig,gf in tmp: tmp[sig].from_L_G_R(self.mat[iNs][iorb],obj[iorb][sig],self.mat[iNs][iorb].conjugate().transpose())
                        tmp *= 1.0/self.Ns
                        symm_obj[jorb] += tmp

                else:

                    if (type(obj[iorb])==DictType):

                        for ii in obj[iorb]:
                            if (l==0):
                                symm_obj[jorb][ii] += obj[iorb][ii]/self.Ns
                            else:
                                symm_obj[jorb][ii] += numpy.dot(numpy.dot(self.mat[iNs][iorb],obj[iorb][ii]),self.mat[iNs][iorb].conjugate().transpose()) / self.Ns

                            

                    else:
                        if (l==0):
                            symm_obj[jorb] += obj[iorb]/self.Ns
                        else:
                            symm_obj[jorb] += numpy.dot(numpy.dot(self.mat[iNs][iorb],obj[iorb]),self.mat[iNs][iorb].conjugate().transpose()) / self.Ns
                        
     
        #for iorb in range(self.N_orbits):
        #    symm_obj[iorb] *= 1.0/self.Ns
    
        return symm_obj
                    
                
    def symmetrise_SO(self,obj):
        """ symmetrisation with SO coupling"""

        if (isinstance(obj[0],GF)):
            symm_obj = [ obj[i].copy() for i in range(len(obj)) ]        # here the result is stored, it is a GF!
            for iorb in range(self.N_orbits): symm_obj[iorb].zero()      # set to zero
        else:
            # if not a GF, we assume it is a matrix (density matrix), has to be complex since self.mat is complex!
            # first index: orbital
            #symm_obj = [ numpy.zeros([self.orbits[iorb][3],self.orbits[iorb][3]],numpy.complex_) for iorb in range(self.N_orbits) ]
            symm_obj = [ copy.deepcopy(obj[i]) for i in range(len(obj)) ]
            
            for iorb in range(self.N_orbits):
                if (type(symm_obj[iorb])==DictType):
                    for ii in symm_obj[iorb]: symm_obj[iorb][ii] *= 0.0
                else:
                    symm_obj[iorb] *= 0.0
                               #print symm_obj[iorb][ii].trace()
        
                
        for iNs in range(self.Ns):

            for iorb in range(self.N_orbits):

                l = self.orbits[iorb][2]         # s, p, d, or f
                dim = self.orbits[iorb][3]
                jorb = self.map[iNs][iorb]
             
                if (isinstance(obj[0],GF)):

                    tmp = obj[iorb].copy()
                    #if self.timeinv[iNs]: tmp <<= tmp.conjugate()
                    if self.timeinv[iNs]: tmp <<= tmp.transpose()

                    for sig,gf in tmp:

                        if (l>0):
                            for sp1 in [0,1]:
                                for sp2 in [0,1] :
                                    tmp[sig][sp1*dim:(sp1+1)*dim,sp2*dim:(sp2+1)*dim].from_L_G_R(self.mat[iNs][iorb],tmp[sig][sp1*dim:(sp1+1)*dim,sp2*dim:(sp2+1)*dim],
                                                                                                 self.mat[iNs][iorb].conjugate().transpose())

                        tmp[sig][0:dim,dim:2*dim] *= (cos(self.phase[iNs])+1j*sin(self.phase[iNs]))
                        tmp[sig][dim:2*dim,0:dim] *= (cos(self.phase[iNs])-1j*sin(self.phase[iNs]))
                  
                    tmp *= 1.0/self.Ns
                    symm_obj[jorb] += tmp

                else:

                    tmp = copy.deepcopy(obj[iorb])

                    if (type(obj[iorb])==DictType):
                        for ii in tmp:
                            if self.timeinv[iNs]: tmp[ii] = tmp[ii].transpose() #tmp[ii].conjugate()

                            if (l>0):
                                for sp1 in [0,1]:
                                    for sp2 in [0,1] :
                                        tmp[ii][sp1*dim:(sp1+1)*dim,sp2*dim:(sp2+1)*dim] = numpy.dot(numpy.dot(self.mat[iNs][iorb],
                                                                                                               tmp[ii][sp1*dim:(sp1+1)*dim,sp2*dim:(sp2+1)*dim]),
                                                                                                     self.mat[iNs][iorb].conjugate().transpose())
			    

                            tmp[ii][0:dim,dim:2*dim] *= (cos(self.phase[iNs])+1j*sin(self.phase[iNs]))
                            tmp[ii][dim:2*dim,0:dim] *= (cos(self.phase[iNs])-1j*sin(self.phase[iNs]))

                            symm_obj[jorb][ii] += tmp[ii]/self.Ns

                    else:
                               
                        if self.timeinv[iNs]: tmp = tmp.transpose()   #tmp.conjugate()

                        if (l>0):
                            for sp1 in [0,1]:
                                for sp2 in [0,1] :
                                    tmp[sp1*dim:(sp1+1)*dim,sp2*dim:(sp2+1)*dim] = numpy.dot(numpy.dot(self.mat[iNs][iorb],
                                                                                                       tmp[sp1*dim:(sp1+1)*dim,sp2*dim:(sp2+1)*dim]),
                                                                                             self.mat[iNs][iorb].conjugate().transpose())

                        tmp[0:dim,dim:2*dim] *= (cos(self.phase[iNs])+1j*sin(self.phase[iNs]))
                        tmp[dim:2*dim,0:dim] *= (cos(self.phase[iNs])-1j*sin(self.phase[iNs]))
                    
                        symm_obj[jorb] += tmp/self.Ns
            
        #for iorb in range(self.N_orbits):
        #    symm_obj[iorb] *= 1.0/self.Ns

    
        return symm_obj
                    
                


    def symmetrise(self,obj):
        """ this should symmetrise the given object
        """

        assert isinstance(obj,list),"obj has to be a list of objects!"
        assert len(obj)==self.N_orbits,"obj has to be a list of the same length as defined in the init"

        if (self.SO==0): 
            symmobj = self.symmetrise_noSO(obj=obj)
        else:
            symmobj = self.symmetrise_SO(obj=obj)
    

        return symmobj

            

        
