
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

from types import *
from pytriqs.dft.symmetry_so import *
import numpy
import pytriqs.base.utility.Dichotomy as Dichotomy
from pytriqs.base.gf_local.GF import GF
from pytriqs.base.gf_local.GFBloc_ImFreq import GFBloc_ImFreq
from pytriqs.base.gf_local.GFBloc_ReFreq import GFBloc_ReFreq
from pytriqs.base.gf_local import GF_Initializers
from pytriqs.solvers.operators import *
from pytriqs.base.archive.HDF_Archive import *
import pytriqs.base.utility.MPI as MPI

from math import cos,sin

import string, pickle


def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
	    yield string.atof(x)


class SumK_LDA_SO:
    """This class provides a general SumK method for combining ab-initio code and pytriqs. 
    At initialization, the things for the DMFT SC loop are read, nothing more!
    
    Several things have to be provided in the LDA data file:
    The LDA Hamiltonian, eps(k)
    The Projection Matrix
    For the k-sum: k-points, weights, and the Rotation Matrix for the symmetry operations

    
    N_Orbitals: Number of orbitals of the embedded G and Sigma = Energy window used in LDA!
                In this implementation, the number of bands can depend on the k point!
    N_Orbitals_corr: Number of correlated orbitals per atom, independent of k
    N_Atoms: Number of atoms with corr. orbitals in the unit cell (index R)
    k_dep_projection: determines if projection matrix is k dependent, 0 means independent
    EnergyUnit: Multiplicative factor to get energies in eV
    N_symm_operations: Number of symmetries used for k sum

    Up to now, no spin dependence!!
    """


    def __init__(self, LDA_file, Symm_file = 'symmetries.dat', mu = 0, UseLDABlocs = False):
        """
        This init differs from SumK_LDA, it contains the modifications for the SO case.
        """
        
        # check the given input:
        assert type(LDA_file)==StringType,"LDA_file must be a filename"; self.LDA_file = LDA_file
       
        self.Chemical_Potential = mu
        
        self.__G_calculated = False                                      # a flag so that we know if we have already done the sum over k
        self.__Proj_Mat_pc_read = False    # Set a flag if Density Projectors to local orbitals are read

        # R is a generator : each R.Next() will return the next number in the file
        R = Read_Fortran_File(LDA_file)
        try:
            self.EnergyUnit = R.next()                         # read the energy convertion factor
            self.Nk = int(R.next())                            # read the number of k points
            self.k_dep_projection = int(R.next())              # read the flag for k-dep. projection matrices
            self.SP = int(R.next())                            # flag for spin-polarised calculation
            self.SO = int(R.next())                            # flag for spin-orbit calculation
            self.charge_below = R.next()                       # total charge below energy window
            self.Density_Required = R.next()                   # total density required, for setting the chemical potential
            self.symm_op = int(R.next())                         # Use symmetry groups for the k-sum?

            # the information on the non-correlated shells is not important here, maybe skip:
            self.N_shells = int(R.next())                      # number of shells (e.g. Fe d, As p, O p) in the unit cell, 
                                                               # corresponds to index R in formulas
            # now read the information about the shells:
            #if (self.SO==0):
            #    self.shells = [ [ int(R.next()) for i in range(4) ] for icrsh in range(self.N_shells) ]    # reads iatom, sort, l, dim
            #else:
            self.shells = [ [ int(R.next()) for i in range(4) ] for icrsh in range(self.N_shells) ]    # reads iatom, sort, l, dim

            self.N_corr_shells = int(R.next())                 # number of corr. shells (e.g. Fe d, Ce f) in the unit cell, 
                                                               # corresponds to index R in formulas
            # now read the information about the shells:
            #if (self.SO==0):
            #    self.corr_shells = [ [ int(R.next()) for i in range(4) ] for icrsh in range(self.N_corr_shells) ]    # reads iatom, sort, l, dim
            #else:
            self.corr_shells = [ [ int(R.next()) for i in range(5) ] for icrsh in range(self.N_corr_shells) ]    # reads iatom, sort, l, dim, SO flag

            self.inequiv_shells(self.corr_shells)              # determine the number of inequivalent correlated shells


            self.use_rotations = int(R.next())
            self.rotmat = [numpy.identity(self.corr_shells[icrsh][3],numpy.complex_) for icrsh in xrange(self.N_corr_shells)]
            if (self.use_rotations):
                # read the matrices
                if (self.corr_shells[icrsh][4]==1): self.rotmat_ph, self.rotmat_timeinv = range(self.N_corr_shells), range(self.N_corr_shells)

                for icrsh in xrange(self.N_corr_shells):
                    for i in xrange(self.corr_shells[icrsh][3]):    # read real part:
                        for j in xrange(self.corr_shells[icrsh][3]):
                            self.rotmat[icrsh][j,i] = R.next()
                    for i in xrange(self.corr_shells[icrsh][3]):    # read imaginary part:
                        for j in xrange(self.corr_shells[icrsh][3]):
                            self.rotmat[icrsh][j,i] += 1j * R.next()

                    if (self.corr_shells[icrsh][4]==1):             # read SO props:
                        self.rotmat_ph[icrsh] = R.next()
                        self.rotmat_timeinv[icrsh] = int(R.next())
                    
                  
            
            # Read here the transformation matrix for complex->cubic:
            self.T = [numpy.zeros([2*self.corr_shells[self.invshellmap[ish]][2]+1,2*self.corr_shells[self.invshellmap[ish]][2]+1],numpy.complex_) 
                      for ish in xrange(self.N_inequiv_corr_shells)]

            for ish in xrange(self.N_inequiv_corr_shells):
                ll = 2*self.corr_shells[self.invshellmap[ish]][2]+1
                for i in xrange(ll):
                    for j in xrange(ll):
                        self.T[ish][j,i] = R.next()
                for i in xrange(ll):
                    for j in xrange(ll):
                        self.T[ish][j,i] += 1j * R.next()
    
            # Spin blocks to be read:
            #if (self.SP==1):
            #    self.Nspinblocs = 2
            #else:
            #    self.Nspinblocs = 1
            self.Nspinblocs = self.SP + 1 - self.SO   # number of spins to read for Norbs and Ham, NOT Projectors
                  

            # GF structure used for the local things in the k sums
            #if (self.SO==0):
            #    self.blocnames = ['up','down']
            #else:
            #    self.blocnames = ['ud']
            #self.NspinblocsGF = len(self.blocnames)
            self.blocnames= [ ['up','down'], ['ud'] ]
            self.NspinblocsGF = [2,1]

            self.names_to_ind = [{}, {}]
            for ibl in range(2):
                for inm in range(self.NspinblocsGF[ibl]): 
                    self.names_to_ind[ibl][self.blocnames[ibl][inm]] = inm * self.SP #(self.Nspinblocs-1)

            #self.GFStruct_corr = [ [ ('up', range(self.corr_shells[i][3])), ('down', range(self.corr_shells[i][3])) ]  for i in xrange(self.N_corr_shells) ]
            self.GFStruct_corr = [ [ (al, range( ((self.corr_shells[i][4]+1) * self.corr_shells[i][3]))) for al in self.blocnames[self.corr_shells[i][4]] ]  
                                   for i in xrange(self.N_corr_shells) ]

            self.GFStruct_Solver = [ [ (al, range( self.corr_shells[self.invshellmap[i]][3]*(1+self.corr_shells[self.invshellmap[i]][4])) ) 
                                       for al in self.blocnames[self.corr_shells[self.invshellmap[i]][4]] ]
                                     for i in xrange(self.N_inequiv_corr_shells) ]
            self.map = [ {} for i in xrange(self.N_inequiv_corr_shells) ]
            self.mapinv = [ {} for i in xrange(self.N_inequiv_corr_shells) ]
            for i in xrange(self.N_inequiv_corr_shells):
                for al in self.blocnames[self.corr_shells[self.invshellmap[i]][4]]:
                    self.map[i][al] = [al for j in range( self.corr_shells[self.invshellmap[i]][3]*(1+self.corr_shells[self.invshellmap[i]][4]) ) ]
                    self.mapinv[i][al] = al
            

            if (self.symm_op==1):
                # Symmetries are used, initialise Symmetry class:
                self.Symm_corr = Symmetry(Symm_file, self.corr_shells, SO = self.SO)
                
            
            
            if self.k_dep_projection == 0:
                # Projection Matrix does not depend on k:
                # i.e., number of bands in the energy window is the same for all k!!!
                # read it once and for all

                no = int(R.next())                             # Orbitals in the energy window of LDA
                self.N_Orbitals = [ [no for isp in range(self.Nspinblocs)] for ik in xrange(self.Nk)]

                # Initialise P, a list of matrices:
                self.Proj_Mat = [ [numpy.zeros([self.corr_shells[icrsh][3], no], numpy.complex_) for icrsh in range (self.N_corr_shells) ]
                                  for isp in range(self.Nspinblocs) ]

                for isp in range(self.Nspinblocs):
                    for icrsh in range(self.N_corr_shells):
                                        
                        for i in xrange(self.corr_shells[icrsh][3]):    # read real part:
                            for j in xrange(no):
                                self.Proj_Mat[isp][icrsh][i,j] = R.next()

                        for i in xrange(self.corr_shells[icrsh][3]):    # read imaginary part:
                            for j in xrange(no):
                                self.Proj_Mat[isp][icrsh][i,j] += 1j * R.next()

            else:
                # read the list of N_Orbitals for all k points
                self.N_Orbitals = [ [0 for isp in range(self.Nspinblocs)] for ik in xrange(self.Nk)]
                for isp in range(self.Nspinblocs):
                    for ik in xrange(self.Nk):
                        self.N_Orbitals[ik][isp] = int(R.next())
                #self.N_Orbitals = [ [int(R.next()) for isp in range(self.Nspinblocs)] for ik in range(self.Nk)]      

                # Initialise P, here a double list of matrices:
                spmap = [[0,1],[0,0]]    # mapping of spin number for Norbitals
                # in SO case, only one projector is constructed:
                self.Proj_Mat = [ [ [numpy.zeros([self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4]), self.N_Orbitals[ik][isp]], numpy.complex_) 
                                     for icrsh in range (self.N_corr_shells)] 
                                    for isp in range(self.Nspinblocs)] 
                                  for ik in range(self.Nk) ]
 
                for isp in range(self.SP+1): # NOTE: in SO case, we have self.SP=1!
                    for ik in xrange(self.Nk):
                        for icrsh in range(self.N_corr_shells):
                            no = self.corr_shells[icrsh][3]

                            for i in xrange(no):    # read real part:
                                for j in xrange(self.N_Orbitals[ik][spmap[self.SO][isp]]):
                                    if (self.SO==0):
                                        self.Proj_Mat[ik][isp][icrsh][i,j] = R.next()
                                    else:
                                        self.Proj_Mat[ik][0][icrsh][isp*no+i,j] = R.next()
                                
                            for i in xrange(no):    # read imaginary part:
                                for j in xrange(self.N_Orbitals[ik][spmap[self.SO][isp]]):
                                    if (self.SO==0):
                                        self.Proj_Mat[ik][isp][icrsh][i,j] += 1j * R.next()
                                    else:
                                        self.Proj_Mat[ik][0][icrsh][isp*no+i,j] += 1j * R.next()

            
                
          
            # now define the arrays for weights and hopping ...
            self.BZ_weights = numpy.ones([self.Nk],numpy.float_)/ float(self.Nk)  # w(k_index),  default normalisation 
            self.Hopping = [ [numpy.zeros([self.N_Orbitals[ik][isp],self.N_Orbitals[ik][isp]],numpy.complex_) 
                              for isp in range(self.Nspinblocs)] for ik in xrange(self.Nk) ]

            # read all the stuff:
            
            # weights in the file
            for ik in xrange(self.Nk) : 
                self.BZ_weights[ik] = R.next()         
                

            # if the sum over spins is in the weights, take it out again!!
            sm = sum(self.BZ_weights)
            self.BZ_weights[:] /= sm 
	    
            # Grab the H
            for isp in range(self.Nspinblocs):
                for ik in xrange(self.Nk) :

                    no = self.N_Orbitals[ik][isp]

                    for i in xrange(no) :
                        for j in xrange(i,no) :
                            self.Hopping[ik][isp][i,j] = R.next() * self.EnergyUnit
		       
                    for i in xrange(no) :
                        for j in xrange(i,no) :
                            self.Hopping[ik][isp][i,j] += 1j * R.next() * self.EnergyUnit
                            if j !=i : self.Hopping[ik][isp][j,i] = self.Hopping[ik][isp][i,j].conjugate()
                            
          
        except StopIteration : # a more explicit error if the file is corrupted.
	    raise "SumK_LDA : reading file HMLT_file failed!"
        
        R.close()

        self.DCenerg = [0.0 for i in xrange(self.N_corr_shells)]

        if (UseLDABlocs):
            dm=self.analyse_BS()
        #else:
            # define a standard GFStruct_Solver and mapping:
            #self.GFStruct_Solver = [ [ ('up', range(self.corr_shells[self.invshellmap[i]][3])), ('down', range(self.corr_shells[self.invshellmap[i]][3])) ]  
            #                           for i in xrange(self.N_inequiv_corr_shells) ]
            #self.map = [ {'down' : ['down' for j in range(self.corr_shells[self.invshellmap[i]][3])], 
            #              'up' : ['up' for j in range(self.corr_shells[self.invshellmap[i]][3])] } for i in xrange(self.N_inequiv_corr_shells) ]
            #self.mapinv = [ {'down' : 'down', 'up' : 'up'} for i in xrange(self.N_inequiv_corr_shells) ]
            

    def save_HDF(self,arxiv):
        """Saves some quantities into an HDF5 arxiv"""

        assert isinstance(arxiv,HDF_Archive),"Give an HDF Archive"

        arxiv['SumK_LDA'] = {'Chemical_Potential' : self.Chemical_Potential,
                             'GFStruct_Solver' : self.GFStruct_Solver,
                             'map' : self.map,
                             'mapinv' : self.mapinv}
        if hasattr(self,"dc_imp"):
            arxiv['SumK_LDA']['dc_imp'] = self.dc_imp
            

    def load_HDF(self,arxiv):
        """Loads some quantities from an HDF5 arxiv"""

        assert isinstance(arxiv,HDF_Archive),"Give an HDF Archive"
        try:
            self.Chemical_Potential = arxiv['SumK_LDA']['Chemical_Potential']
            self.GFStruct_Solver = arxiv['SumK_LDA']['GFStruct_Solver']
            self.map = [ {} for i in xrange(self.N_inequiv_corr_shells) ]
            self.mapinv = [ {} for i in xrange(self.N_inequiv_corr_shells) ]
            for ii in range(self.N_inequiv_corr_shells):
                for nm in self.blocnames[self.corr_shells[ii][4]]:
                    self.map[ii][nm] = arxiv['SumK_LDA']['map'][ii][nm]
                for nm,k in self.GFStruct_Solver[ii]:
                    self.mapinv[ii][nm] = arxiv['SumK_LDA']['mapinv'][ii][nm]
            #self.mapinv = arxiv['SumK_LDA']['mapinv']
            if 'dc_imp' in arxiv['SumK_LDA']:
                self.__initDC()
                for ii in range(self.N_inequiv_corr_shells):
                    for nm in self.dc_imp[ii]:
                        self.dc_imp[ii][nm] = arxiv['SumK_LDA']['dc_imp'][ii][nm]
                #self.dc_imp = arxiv['SumK_LDA']['dc_imp']
            return True
        except:
            MPI.report("Loading failed, starting from scratch...")
            return False
   

     
    def save(self,Filename):
        """Saves some quantities into a file"""
        
        if not MPI.IS_MASTER_NODE(): return   # Do nothing if not master
        assert type(Filename)==StringType,"Filename must be a filename"
        
        f=open(Filename, 'w')
        
        # use pickle module:
        P=pickle.Pickler(f)
        P.dump(self.Chemical_Potential)
        P.dump(self.GFStruct_Solver)
        P.dump(self.map)
        P.dump(self.mapinv)
        if hasattr(self,"dc_imp"):
            fl = 1
            P.dump(fl)
            P.dump(self.dc_imp)
        else:
            fl=0
            P.dump(fl)

        f.close()

       
    def load(self,Filename):
        """Reads some quantities from a file"""
        import os.path

        assert type(Filename)==StringType,"Filename must be a filename"
        if not(os.path.exists(Filename)) :    #raise IOError, "File %s does not exists"%Filename
            # no file found:
            MPI.report("Loading failed, starting from scratch...")
            #self.analyse_BS()
            return False
        else:
            # use pickle module:
            f=open(Filename,'r')
            P=pickle.Unpickler(f)
            self.Chemical_Potential = P.load()
            self.GFStruct_Solver = P.load()
            self.map = P.load()
            self.mapinv = P.load()
            fl = P.load()
            if (fl==1):
                self.dc_imp = P.load()

            return True


    def load_old(self,Filename):
        """The old version of loading, without the block structure!
           It is still included for compatibility."""

        assert type(Filename)==StringType,"Filename must be a filename"

        R=Read_Fortran_File(Filename)
        self.Chemical_Potential = R.next()
        
        if (int(R.next())):
            # construct dc_imp:
            self.dc_imp = [ {} for i in xrange(self.N_corr_shells)]
            for i in xrange(self.N_corr_shells):
                l = self.corr_shells[i][3]
                for j in xrange(len(self.GFStruct_corr[i])):
                    self.dc_imp[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.identity(l,numpy.float_)
            # read from file:
            for ish in xrange(self.N_corr_shells):
                for bl,indices in self.GFStruct_corr[ish]:
                    for ind in indices:
                        self.dc_imp[ish][bl][ind,ind] = R.next()

        # with the standard version of the python scripts, no LDABloc analysis is done if previous
        # results are loaded, so do it now
        dm=self.analyse_BS()


       


    def downfold(self,ik,icrsh,sig,gf_to_downfold,gf_inp):
        """Downfolding a block of the Greens function"""
        
        gf_downfolded = gf_inp.copy()

        #if (self.corr_shells[icrsh][4]==0):
        #    # no SO coupling:
        #    isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        #    gf_downfolded.from_L_G_R(self.Proj_Mat[ik][isp][icrsh],gf_to_downfold,self.Proj_Mat[ik][isp][icrsh].conjugate().transpose())  # downfolding G
        #else:
        #    Norb = self.corr_shells[icrsh][3]
        #    for sp1 in [0,1]:
        #        for sp2 in [0,1]:
        #            gf_downfolded[sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb].from_L_G_R(self.Proj_Mat[ik][sp1][icrsh],
        #                                                                                  gf_to_downfold,self.Proj_Mat[ik][sp2][icrsh].conjugate().transpose())

        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        gf_downfolded.from_L_G_R(self.Proj_Mat[ik][isp][icrsh],gf_to_downfold,self.Proj_Mat[ik][isp][icrsh].conjugate().transpose())  # downfolding G

        return gf_downfolded
        

    def upfold(self,ik,icrsh,sig,gf_to_upfold,gf_inp):
        """Upfolding a block of the Greens function"""

        gf_upfolded = gf_inp.copy()
        #gf_upfolded.zero()

        #if (self.corr_shells[icrsh][4]==0):
            # no SO coupling
        #    isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        #    gf_upfolded.from_L_G_R(self.Proj_Mat[ik][isp][icrsh].conjugate().transpose(),gf_to_upfold,self.Proj_Mat[ik][isp][icrsh]) 
        #else:
        #    Norb = self.corr_shells[icrsh][3]
        #    gftmp = gf_upfolded.copy()
        #    for sp1 in [0,1]:
        #        for sp2 in [0,1]:
        #            gftmp.from_L_G_R(self.Proj_Mat[ik][sp1][icrsh].conjugate().transpose(),
        #                             gf_to_upfold[sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb],self.Proj_Mat[ik][sp2][icrsh])
        #            gf_upfolded += gftmp

        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        gf_upfolded.from_L_G_R(self.Proj_Mat[ik][isp][icrsh].conjugate().transpose(),gf_to_upfold,self.Proj_Mat[ik][isp][icrsh]) 

        return gf_upfolded


    def rotloc(self,icrsh,gf_to_rotate,direction):
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

        assert ((direction=='toLocal')or(direction=='toGlobal')),"Give direction 'toLocal' or 'toGlobal' in rotloc!"

        gf_rotated = gf_to_rotate.copy()
        if (direction=='toLocal'):
            if (self.corr_shells[icrsh][4]==0):
                gf_rotated.from_L_G_R(self.rotmat[icrsh].transpose(),gf_to_rotate,self.rotmat[icrsh].conjugate())
            else:
                if (self.rotmat_timeinv[icrsh]==1): gf_rotated <<= gf_to_rotate.transpose()
                Norb = self.corr_shells[icrsh][3]
                for sp1 in [0,1]: 
                    for sp2 in [0,1]:
                        gf_rotated[sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb].from_L_G_R(self.rotmat[icrsh].transpose(),gf_rotated[sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb],self.rotmat[icrsh].conjugate())
                
                gf_rotated[0:Norb,Norb:2*Norb] *= (cos(self.rotmat_ph[icrsh])+1j*sin(self.rotmat_ph[icrsh]))
                gf_rotated[Norb:2*Norb,0:Norb] *= (cos(self.rotmat_ph[icrsh])-1j*sin(self.rotmat_ph[icrsh]))

        elif (direction=='toGlobal'):
            if (self.corr_shells[icrsh][4]==0):
                gf_rotated.from_L_G_R(self.rotmat[icrsh].conjugate(),gf_to_rotate,self.rotmat[icrsh].transpose())
            else:
                if (self.rotmat_timeinv[icrsh]==1): gf_rotated <<= gf_to_rotate.transpose()
                Norb = self.corr_shells[icrsh][3]
                for sp1 in [0,1]: 
                    for sp2 in [0,1]:
                        gf_rotated[sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb].from_L_G_R(self.rotmat[icrsh].conjugate(),gf_rotated[sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb],self.rotmat[icrsh].transpose())
                
                gf_rotated[0:Norb,Norb:2*Norb] *= (cos(self.rotmat_ph[icrsh])-1j*sin(self.rotmat_ph[icrsh]))
                gf_rotated[Norb:2*Norb,0:Norb] *= (cos(self.rotmat_ph[icrsh])+1j*sin(self.rotmat_ph[icrsh]))
    
        return gf_rotated
                    


    


    def check_projectors(self):

        densmat = [numpy.zeros([self.corr_shells[ish][3]*(self.corr_shells[ish][4]+1),self.corr_shells[ish][3]*(self.corr_shells[ish][4]+1)],numpy.complex_) 
                   for ish in range(self.N_corr_shells)]
        
        for ik in range(self.Nk):
        
            for ish in range(self.N_corr_shells):
                Norb = self.corr_shells[ish][3]
                #for sp1 in [0,1]:
                #    for sp2 in [0,1]:
                densmat[ish][:,:] += numpy.dot(self.Proj_Mat[ik][0][ish],self.Proj_Mat[ik][0][ish].transpose().conjugate()) * self.BZ_weights[ik]

        if (self.symm_op!=0): densmat = self.Symm_corr.symmetrise(densmat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                if (self.corr_shells[icrsh][4]==0):
                    densmat[icrsh] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),densmat[icrsh]) , 
                                                self.rotmat[icrsh].conjugate() )
                 
                else:
                    if (self.rotmat_timeinv[icrsh]==1): densmat[icrsh] = densmat[icrsh].conjugate()
                    Norb = self.corr_shells[icrsh][3]
                    for sp1 in [0,1]: 
                        for sp2 in [0,1]:
                            densmat[icrsh][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),densmat[icrsh][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb]),self.rotmat[icrsh].conjugate())
                    densmat[icrsh][0:Norb,Norb:2*Norb] *= (cos(self.rotmat_ph[icrsh])+1j*sin(self.rotmat_ph[icrsh]))
                    densmat[icrsh][Norb:2*Norb,0:Norb] *= (cos(self.rotmat_ph[icrsh])-1j*sin(self.rotmat_ph[icrsh]))


         
        return densmat



    def simplepointdensmat(self):


        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        MMat = [numpy.zeros( [self.N_Orbitals[0][ntoi[bl]],self.N_Orbitals[0][ntoi[bl]]], numpy.complex_) for bl in bln] 

        #densmat = [ numpy.zeros([self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4]),self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4])], 
        #                         numpy.complex_) for icrsh in xrange(self.N_corr_shells) ]  

        densmat = [ {} for icrsh in xrange(self.N_corr_shells)]
        for icrsh in xrange(self.N_corr_shells):
            for bl in self.blocnames[self.corr_shells[icrsh][4]]:
                densmat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4]),
                                                  self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4])], numpy.complex_)


        ikarray=numpy.array(range(self.Nk))
          
        for ik in MPI.slice_array(ikarray):
            
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==len(MMat[ib]) 
                                   for ib in range(self.NspinblocsGF[self.SO]) ] )
               
            if (not unchangedsize):
                MMat = [numpy.zeros( [self.N_Orbitals[ik][ntoi[bl]],self.N_Orbitals[ik][ntoi[bl]]], numpy.complex_) for bl in bln] 

            for bl in bln:
                ind = ntoi[bl]
                for inu in range(self.N_Orbitals[ik][ind]):
                    if (self.Hopping[ik][ind][inu,inu] < 0.0): 
                        MMat[ind][inu,inu] = 1.0
                    else:
                        MMat[ind][inu,inu] = 0.0 


            for icrsh in range(self.N_corr_shells):
                for bn in self.blocnames[self.corr_shells[icrsh][4]]:
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    #print ik, bn, isp
                    densmat[icrsh][bn] += self.BZ_weights[ik] * numpy.dot( numpy.dot(self.Proj_Mat[ik][isp][icrsh],MMat[isp]) , 
                                                                           self.Proj_Mat[ik][isp][icrsh].transpose().conjugate() )

        # get data from nodes:
        for icrsh in range(self.N_corr_shells):
            for sig in densmat[icrsh]:
                densmat[icrsh][sig] = MPI.all_reduce(MPI.world,densmat[icrsh][sig],lambda x,y : x+y)
        MPI.barrier()
                    
        if (self.symm_op!=0): densmat = self.Symm_corr.symmetrise(densmat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for bn in densmat[icrsh]:

                    if (self.corr_shells[icrsh][4]==0):
                        densmat[icrsh][bn] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),densmat[icrsh][bn]) , 
                                                        self.rotmat[icrsh].conjugate() )
                 
                    else:
                        if (self.rotmat_timeinv[icrsh]==1): densmat[icrsh][bn] = densmat[icrsh][bn].conjugate()
                        Norb = self.corr_shells[icrsh][3]
                        for sp1 in [0,1]: 
                            for sp2 in [0,1]:
                                densmat[icrsh][bn][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),densmat[icrsh][bn][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb]),self.rotmat[icrsh].conjugate())
                            densmat[icrsh][bn][0:Norb,Norb:2*Norb] *= (cos(self.rotmat_ph[icrsh])+1j*sin(self.rotmat_ph[icrsh]))
                            densmat[icrsh][bn][Norb:2*Norb,0:Norb] *= (cos(self.rotmat_ph[icrsh])-1j*sin(self.rotmat_ph[icrsh]))

       
        return densmat


    def density_gf(self,Beta = 40):
        """Calculates the density without setting up Gloc. It is useful for Hubbard I, and very fast.""" 

        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        densmat = [ {} for icrsh in xrange(self.N_corr_shells)]
        for icrsh in xrange(self.N_corr_shells):
            for bl in self.blocnames[self.corr_shells[icrsh][4]]:
                densmat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4]),
                                                  self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4])], numpy.complex_)

        # initialisation:
        BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
        GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
        a_list = [a for a,al in GFStruct]   
        glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta = Beta) for a,al in GFStruct]  
        Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
        Gupf.zero()
        mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln] 
        for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential


        ikarray=numpy.array(range(self.Nk))
        #ikarray=numpy.array([0])

        for ik in MPI.slice_array(ikarray):
            
            GFsize = [ gf.N1 for sig,gf in Gupf]  
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                   for ib in range(self.NspinblocsGF[self.SO]) ] )
               
                #if (self.N_Orbitals[ik]!=GFsize):
            if (not unchangedsize):
                BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                a_list = [a for a,al in GFStruct]                                 
                glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta = Beta) for a,al in GFStruct]    
                Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
                Gupf.zero()
                mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential
                  

            Gupf <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
            M = copy.deepcopy(mupat)
            for ibl in range(self.NspinblocsGF[self.SO]): 
                ind = ntoi[bln[ibl]]
                M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
            Gupf -= M

            Gupf.invert()
            Gupf *= self.BZ_weights[ik]

            MMat = range(len(bln))
            dm = Gupf.density()

            

            for bl in bln:
                ind = ntoi[bl]
                MMat[ind] = dm[bl]

            
            for icrsh in range(self.N_corr_shells):
                for bn in self.blocnames[self.corr_shells[icrsh][4]]:
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    #print ik, bn, isp
                    densmat[icrsh][bn] += numpy.dot( numpy.dot(self.Proj_Mat[ik][isp][icrsh],MMat[isp]),self.Proj_Mat[ik][isp][icrsh].transpose().conjugate() )


        if (self.symm_op!=0): densmat = self.Symm_corr.symmetrise(densmat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for bn in densmat[icrsh]:
                    if (self.corr_shells[icrsh][4]==0):
                        densmat[icrsh][bn] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),densmat[icrsh][bn]) , 
                                                    self.rotmat[icrsh].conjugate() )
                 
                    else:
                        if (self.rotmat_timeinv[icrsh]==1): densmat[icrsh][bn] = densmat[icrsh][bn].conjugate()
                        Norb = self.corr_shells[icrsh][3]
                        for sp1 in [0,1]: 
                            for sp2 in [0,1]:
                                densmat[icrsh][bn][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),densmat[icrsh][bn][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb]),self.rotmat[icrsh].conjugate())
                        densmat[icrsh][bn][0:Norb,Norb:2*Norb] *= (cos(self.rotmat_ph[icrsh])+1j*sin(self.rotmat_ph[icrsh]))
                        densmat[icrsh][bn][Norb:2*Norb,0:Norb] *= (cos(self.rotmat_ph[icrsh])-1j*sin(self.rotmat_ph[icrsh]))

        
        dm = [ densmat[self.invshellmap[ish]] for ish in range(self.N_inequiv_corr_shells) ]

        return dm

        
           
    def analyse_BS_from_GF(self,Beta = 40, threshold = 0.00000001, includeshells = None):
        """Analyses the LDA Green function and gives the optimal block sizes for the CTQMC Solver.
           It is done at a given Beta, which is not important in this case. It can differ from the Beta used later for the calculations.
           includeshells can be a list of the inequivalent shells to be included in this analysis, excluded shells will use the standard bloc structure."""

        if not (includeshells is None):
            assert (len(includeshells)<=self.N_inequiv_corr_shells), "Too long list (includeshells) in analyse_BS"
            assert ((max(includeshells)<=self.N_inequiv_corr_shells)and(min(includeshells)>=0)), "Wrong item in list includeshells!"
        
                    
        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        # initialisation:
        BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
        GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
        a_list = [a for a,al in GFStruct]   
        glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta = Beta) for a,al in GFStruct]  
        Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
        Gupf.zero()
        mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln] 
        for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential
        Gloc = [ GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = Gupf.mesh)) for a,al in self.GFStruct_corr[icrsh] ],
                    Copy = False) for icrsh in xrange(self.N_corr_shells) ]   
        for icrsh in xrange(self.N_corr_shells): Gloc[icrsh].zero()                        # initialize to zero


        #for ik in xrange(self.Nk):
        ikarray=numpy.array(range(self.Nk))
        #ikarray=numpy.array([100])

        for ik in MPI.slice_array(ikarray):
            print ik

            GFsize = [ gf.N1 for sig,gf in Gupf]  
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                   for ib in range(self.NspinblocsGF[self.SO]) ] )
               
                #if (self.N_Orbitals[ik]!=GFsize):
            if (not unchangedsize):
                BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                a_list = [a for a,al in GFStruct]                                 
                glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta = Beta) for a,al in GFStruct]    
                Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
                Gupf.zero()
                mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential
                  

            Gupf <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
            M = copy.deepcopy(mupat)
            for ibl in range(self.NspinblocsGF[self.SO]): 
                ind = ntoi[bln[ibl]]
                M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
            Gupf -= M

            Gupf.invert()
            Gupf *= self.BZ_weights[ik]
#            Gfdat = Gupf['ud']._data.array[:,:,0]
#            Gupf.save('gupf')
                
            for icrsh in xrange(self.N_corr_shells):
                tmp = Gloc[icrsh].copy()      # init temporary storage
                for sig,gf in tmp: tmp[sig] <<= self.downfold(ik,icrsh,sig,Gupf[sig],gf)
                    
                Gloc[icrsh] += tmp

            
            #collect data from MPI:
        for icrsh in xrange(self.N_corr_shells):
            Gloc[icrsh] <<= MPI.all_reduce(MPI.world,Gloc[icrsh],lambda x,y : x+y)
        MPI.barrier()

        # output for tests:
        #Gloc[0].save('Glocbfsymm1')
        #Gloc[1].save('Glocbfsymm2')
        #densmat2 = [ Gloc[self.invshellmap[ish]].density() for ish in xrange(self.N_inequiv_corr_shells) ]

        #Gfdat = Gloc[0]['ud']._data.array[:,:,0]

        # Gloc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:   
        if (self.symm_op!=0): Gloc = self.Symm_corr.symmetrise(Gloc)
        
        # output for tests:
        #Gloc[0].save('Glocafsymm1')
        #Gloc[1].save('Glocafsymm2')
        #densmat1 = [ Gloc[self.invshellmap[ish]].density() for ish in xrange(self.N_inequiv_corr_shells) ]


        # Gloc is rotated to the local coordinate system:

        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in Gloc[icrsh]: Gloc[icrsh][sig] <<= self.rotloc(icrsh,gf,direction='toLocal')
        
        # output for tests:
        #Gloc[0].save('Gloc2afrot1')
        #Gloc[1].save('Gloc2afrot2')

        # now collect the density matrices for inequivalent shells:
        densmat = [ Gloc[self.invshellmap[ish]].density() for ish in xrange(self.N_inequiv_corr_shells) ]
        #self.GFStruct_Solver, self.map  = [], [{} for ish in xrange(self.N_inequiv_corr_shells)]
        #self.mapinv = [{} for ish in xrange(self.N_inequiv_corr_shells)]

        #for ish in xrange(self.N_inequiv_corr_shells):
        if includeshells is None: includeshells=range(self.N_inequiv_corr_shells)
        for ish in includeshells:

            #self.GFStruct_Solver.append([])
            self.GFStruct_Solver[ish] = []

            a_list = [a for a,al in self.GFStruct_corr[self.invshellmap[ish]] ]
            for a in a_list:
                
                dm = densmat[ish][a]            
                dmbool = (abs(dm) > threshold)          # gives an index list of entries larger that threshold

                offdiag = []
                for i in xrange(len(dmbool)):
                    for j in xrange(i,len(dmbool)):
                        if ((dmbool[i,j])&(i!=j)): offdiag.append([i,j])

                NBlocs = len(dmbool)
                blocs = [ [i] for i in range(NBlocs) ]

                for i in range(len(offdiag)):
                    if (offdiag[i][0]!=offdiag[i][1]):
                        for j in range(len(blocs[offdiag[i][1]])): blocs[offdiag[i][0]].append(blocs[offdiag[i][1]][j])
                        del blocs[offdiag[i][1]]
                        for j in range(i+1,len(offdiag)):
                            if (offdiag[j][0]==offdiag[i][1]): offdiag[j][0]=offdiag[i][0]
                            if (offdiag[j][1]==offdiag[i][1]): offdiag[j][1]=offdiag[i][0]
                            if (offdiag[j][0]>offdiag[i][1]): offdiag[j][0] -= 1
                            if (offdiag[j][1]>offdiag[i][1]): offdiag[j][1] -= 1
                            offdiag[j].sort()
                        NBlocs-=1

                for i in range(NBlocs):
                    blocs[i].sort()
                    self.GFStruct_Solver[ish].append( ('%s%s'%(a,i),blocs[i]) )
                   
                               
                # map is the mapping of the blocs from the SK blocs to the CTQMC blocs:
                self.map[ish][a] = range(len(dmbool))
                for ibl in range(NBlocs):
                    for j in range(len(blocs[ibl])):
                        self.map[ish][a][blocs[ibl][j]] = '%s%s'%(a,ibl)
                        self.mapinv[ish]['%s%s'%(a,ibl)] = a


        # return the density matrix for the double counting!
        return densmat


    def analyse_BS(self, threshold = 0.00000001, includeshells = None):
        """ Determines the Block structure from the simple point integration"""

        dm = self.simplepointdensmat()
        
        densmat = [dm[self.invshellmap[ish]] for ish in xrange(self.N_inequiv_corr_shells) ]

        if includeshells is None: includeshells=range(self.N_inequiv_corr_shells)
        for ish in includeshells:

            #self.GFStruct_Solver.append([])
            self.GFStruct_Solver[ish] = []

            a_list = [a for a,al in self.GFStruct_corr[self.invshellmap[ish]] ]
            for a in a_list:
                
                dm = densmat[ish][a]            
                dmbool = (abs(dm) > threshold)          # gives an index list of entries larger that threshold

                offdiag = []
                for i in xrange(len(dmbool)):
                    for j in xrange(i,len(dmbool)):
                        if ((dmbool[i,j])&(i!=j)): offdiag.append([i,j])

                NBlocs = len(dmbool)
                blocs = [ [i] for i in range(NBlocs) ]

                for i in range(len(offdiag)):
                    if (offdiag[i][0]!=offdiag[i][1]):
                        for j in range(len(blocs[offdiag[i][1]])): blocs[offdiag[i][0]].append(blocs[offdiag[i][1]][j])
                        del blocs[offdiag[i][1]]
                        for j in range(i+1,len(offdiag)):
                            if (offdiag[j][0]==offdiag[i][1]): offdiag[j][0]=offdiag[i][0]
                            if (offdiag[j][1]==offdiag[i][1]): offdiag[j][1]=offdiag[i][0]
                            if (offdiag[j][0]>offdiag[i][1]): offdiag[j][0] -= 1
                            if (offdiag[j][1]>offdiag[i][1]): offdiag[j][1] -= 1
                            offdiag[j].sort()
                        NBlocs-=1

                for i in range(NBlocs):
                    blocs[i].sort()
                    self.GFStruct_Solver[ish].append( ('%s%s'%(a,i),blocs[i]) )
                   
                               
                # map is the mapping of the blocs from the SK blocs to the CTQMC blocs:
                self.map[ish][a] = range(len(dmbool))
                for ibl in range(NBlocs):
                    for j in range(len(blocs[ibl])):
                        self.map[ish][a][blocs[ibl][j]] = '%s%s'%(a,ibl)
                        self.mapinv[ish]['%s%s'%(a,ibl)] = a

        return densmat
        
    

    def eff_atomic_levels(self):
        """Calculates the effective atomic levels needed as input for the Hubbard I Solver."""

        # define matrices for inequivalent shells:
        eff_atlevels = [ {} for ish in range(self.N_inequiv_corr_shells) ]
        for ish in range(self.N_inequiv_corr_shells):
            for bn in self.blocnames[self.corr_shells[self.invshellmap[ish]][4]]:
                eff_atlevels[ish][bn] = numpy.identity(self.corr_shells[self.invshellmap[ish]][3]*(1+self.corr_shells[self.invshellmap[ish]][4]), numpy.complex_)
 
        # Chemical Potential:
        for ish in xrange(self.N_inequiv_corr_shells): 
            for ii in eff_atlevels[ish]: eff_atlevels[ish][ii] *= -self.Chemical_Potential
        
        # double counting term:
        if hasattr(self,"dc_imp"):
            for ish in xrange(self.N_inequiv_corr_shells): 
                for ii in eff_atlevels[ish]:
                    eff_atlevels[ish][ii] -= self.dc_imp[self.invshellmap[ish]][ii]

        # sum over k:
        if not hasattr(self,"Hsumk"):
            # calculate the sum over k. Does not depend on mu, so do it only once:
            self.Hsumk = [ {} for ish in range(self.N_corr_shells) ]
            for icrsh in range(self.N_corr_shells):
                for bn in self.blocnames[self.corr_shells[icrsh][4]]: 
                    dim = self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4])
                    self.Hsumk[icrsh][bn] = numpy.zeros([dim,dim],numpy.complex_)
          
            for icrsh in range(self.N_corr_shells):
                for bn in self.blocnames[self.corr_shells[icrsh][4]]:
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]

                    for ik in xrange(self.Nk):
                        self.Hsumk[icrsh][bn] += self.BZ_weights[ik] * numpy.dot( numpy.dot(self.Proj_Mat[ik][isp][icrsh],self.Hopping[ik][isp]) , 
                                                                                  self.Proj_Mat[ik][isp][icrsh].conjugate().transpose() )

            # symmetrisation:
            if (self.symm_op!=0): self.Hsumk = self.Symm_corr.symmetrise(self.Hsumk)

            # Rotate to local coordinate system:
            if (self.use_rotations):
                for icrsh in xrange(self.N_corr_shells):
                    for bn in self.Hsumk[icrsh]:

                        if (self.corr_shells[icrsh][4]==0):
                            self.Hsumk[icrsh][bn] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),self.Hsumk[icrsh][bn]) , 
                                                               self.rotmat[icrsh].conjugate() )
                 
                        else:
                            if (self.rotmat_timeinv[icrsh]==1): self.Hsumk[icrsh][bn] = self.Hsumk[icrsh][bn].conjugate()
                            Norb = self.corr_shells[icrsh][3]
                            for sp1 in [0,1]: 
                                for sp2 in [0,1]:
                                    self.Hsumk[icrsh][bn][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb] = numpy.dot( numpy.dot(self.rotmat[icrsh].transpose(),self.Hsumk[icrsh][bn][sp1*Norb:(sp1+1)*Norb,sp2*Norb:(sp2+1)*Norb]),
                                                                                                                    self.rotmat[icrsh].conjugate())
                            self.Hsumk[icrsh][bn][0:Norb,Norb:2*Norb] *= (cos(self.rotmat_ph[icrsh])+1j*sin(self.rotmat_ph[icrsh]))
                            self.Hsumk[icrsh][bn][Norb:2*Norb,0:Norb] *= (cos(self.rotmat_ph[icrsh])-1j*sin(self.rotmat_ph[icrsh]))

           
        # add to matrix:
        for ish in xrange(self.N_inequiv_corr_shells): 
            for bn in eff_atlevels[ish]:
                eff_atlevels[ish][bn] += self.Hsumk[self.invshellmap[ish]][bn]


        return eff_atlevels



    def __initDC(self):
        
        # construct the density matrix dm_imp and double counting arrays
        #self.dm_imp = [ {} for i in xrange(self.N_corr_shells)]
        self.dc_imp = [ {} for i in xrange(self.N_corr_shells)]
        for i in xrange(self.N_corr_shells):
            l = (self.corr_shells[i][4]+1) * self.corr_shells[i][3]
            for j in xrange(len(self.GFStruct_corr[i])):
                #self.dm_imp[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.zeros([l,l],numpy.float_)
                self.dc_imp[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.identity(l,numpy.float_)
        self.DCenerg = [0.0 for i in xrange(self.N_corr_shells)]


    def SetDC_Lichtenstein(self,sigimp):
        """Sets a double counting term according to Lichtenstein et al. PRL2001"""

        assert isinstance(sigimp,list), "sigimp has to be a list of impurity self energies for the correlated shells, even if it is of length 1!"
        assert len(sigimp)==self.N_inequiv_corr_shells, "give exactly one Sigma for each inequivalent corr. shell!"

        self.__initDC()

 
        # transform the CTQMC blocks to the full matrix:
        for icrsh in xrange(self.N_corr_shells):
            s = self.shellmap[icrsh]    # s is the index of the inequivalent shell corresponding to icrsh
            for ibl in range(len(self.GFStruct_Solver[s])):
                for i in range(len(self.GFStruct_Solver[s][ibl][1])):
                    for j in range(len(self.GFStruct_Solver[s][ibl][1])):
                        bl   = self.GFStruct_Solver[s][ibl][0]
                        ind1 = self.GFStruct_Solver[s][ibl][1][i]
                        ind2 = self.GFStruct_Solver[s][ibl][1][j]
                        self.dm_imp[icrsh][self.mapinv[s][bl]][ind1,ind2] = sigimp[s][bl]._data.array[i,j,0].real 
                        # self energy at smallest matsubara, could be done better


        for icrsh in xrange(self.N_corr_shells):
            # trace:
            Sigtr = 0.0
            a_list = [a for a,al in self.GFStruct_corr[icrsh]]
            for bl in a_list: Sigtr += self.dm_imp[icrsh][bl].trace()
            for bl in a_list: self.dc_imp[icrsh][bl] *= (Sigtr / (self.corr_shells[icrsh][3] * 2.0))

            #self.dc_imp[icrsh]['up'][:,:] += self.dc_imp[icrsh]['down'][:,:]
            #self.dc_imp[icrsh]['up'][:,:] /= 2.0
            #self.dc_imp[icrsh]['down'][:,:] = self.dc_imp[icrsh]['up'][:,:]

            

    def SetDoubleCounting(self,densmat,U_interact,J_Hund,orb=0,useDCformula=0,useval=None):
        """Sets the double counting term for inequiv orbital orb
           useDCformula=0: LDA+U FLL double counting, useDCformula=1: Held's formula. 
           useDCformula=2: AMF
           Be sure that you use the correct interaction Hamiltonian!"""
        

        if (not hasattr(self,"dc_imp")): self.__initDC()
                    
                
        dm = [ {} for i in xrange(self.N_corr_shells)]
        for i in xrange(self.N_corr_shells):
            l = self.corr_shells[i][3]*(1+self.corr_shells[i][4])
            for j in xrange(len(self.GFStruct_corr[i])):
                dm[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.zeros([l,l],numpy.float_)
        

        for icrsh in xrange(self.N_corr_shells):

            iorb = self.shellmap[icrsh]    # iorb is the index of the inequivalent shell corresponding to icrsh

            if (iorb==orb):
                # do this orbital

                l = self.corr_shells[icrsh][3]*(1+self.corr_shells[icrsh][4])
                for j in xrange(len(self.GFStruct_corr[icrsh])):
                    self.dc_imp[icrsh]['%s'%self.GFStruct_corr[icrsh][j][0]] = numpy.identity(l,numpy.float_)


                # transform the CTQMC blocks to the full matrix:
                for ibl in range(len(self.GFStruct_Solver[iorb])):
                    for i in range(len(self.GFStruct_Solver[iorb][ibl][1])):
                        for j in range(len(self.GFStruct_Solver[iorb][ibl][1])):
                            bl   = self.GFStruct_Solver[iorb][ibl][0]
                            ind1 = self.GFStruct_Solver[iorb][ibl][1][i]
                            ind2 = self.GFStruct_Solver[iorb][ibl][1][j]
                            dm[icrsh][self.mapinv[iorb][bl]][ind1,ind2] = densmat[bl][i,j]

                M = self.corr_shells[icrsh][3]
                Ncr = {}
                Ncrtot = 0.0
                a_list = [a for a,al in self.GFStruct_corr[icrsh]]
                for bl in a_list:
                    Ncr[bl] = dm[icrsh][bl].trace()
                    Ncrtot += Ncr[bl]
                

                if (useval is None):
                              
                    if (useDCformula==0):
                        self.DCenerg[icrsh] = U_interact / 2.0 * Ncrtot * (Ncrtot-1.0)
                        for bl in a_list:
                            Uav = U_interact*(Ncrtot-0.5) - J_Hund*(Ncr[bl] - 0.5)
                            self.dc_imp[icrsh][bl] *= Uav                              
                            self.DCenerg[icrsh]  -= J_Hund / 2.0 * (Ncr[bl]) * (Ncr[bl]-1.0)
                            MPI.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())
                    elif (useDCformula==1):
                        self.DCenerg[icrsh] = (U_interact + J_Hund * (2.0-(M-1)) / (2*M-1)  ) / 2.0 * Ncrtot * (Ncrtot-1.0)
                        for bl in a_list:
                            # Held's formula, with U_interact the interorbital onsite interaction
                            Uav = (U_interact + J_Hund * (2.0-(M-1)) / (2*M-1)  ) * (Ncrtot-0.5)
                            self.dc_imp[icrsh][bl] *= Uav 
                            MPI.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())
                    elif (useDCformula==2):
                        self.DCenerg[icrsh] = 0.5 * U_interact * Ncrtot * Ncrtot
                        for bl in a_list:
                            # AMF
                            Uav = U_interact*(Ncrtot - Ncr[bl]/M) - J_Hund * (Ncr[bl] - Ncr[bl]/M)
                            self.dc_imp[icrsh][bl] *= Uav
                            self.DCenerg[icrsh] -= (U_interact + (M-1)*J_Hund)/M * 0.5 * Ncr[bl] * Ncr[bl]
                            MPI.report("DC for shell %(icrsh)i and block %(bl)s = %(Uav)f"%locals())
                    
                    # output:
                    MPI.report("DC energy for shell %s = %s"%(icrsh,self.DCenerg[icrsh]))

                else:    
            
                    a_list = [a for a,al in self.GFStruct_corr[icrsh]]
                    for bl in a_list:
                        self.dc_imp[icrsh][bl] *= useval
                    
                    self.DCenerg[icrsh] = useval * Ncrtot

                    # output:
                    MPI.report("DC for shell %(icrsh)i = %(useval)f"%locals())
                    MPI.report("DC energy = %s"%self.DCenerg[icrsh])

        


    def find_DC(self,orb,guess,densmat,densreq=None,precision=0.01):
        """searches for DC in order to fulfill charge neutrality.
           If densreq is given, then DC is set such that the LOCAL charge of orbital
           orb coincides with densreq."""
        
        mu = self.Chemical_Potential
        
        def F(dc):
            self.SetDoubleCounting(densmat=densmat,U_interact=0,J_Hund=0,orb=orb,useval=dc)
            if (densreq is None):
                return self.total_density(mu=mu)
            else:
                return self.extract_Gloc()[orb].total_density()
          

        if (densreq is None):
            Dens_rel = self.Density_Required - self.charge_below
        else:
            Dens_rel = densreq
        
        dcnew = Dichotomy.Dichotomy(Function = F,
                                    xinit = guess, yvalue = Dens_rel,
                                    Precision_on_y = precision, Delta_x=0.5,
                                    MaxNbreLoop = 100, xname="Double-Counting", yname= "Total Density",
                                    verbosity = 3)[0]

        return dcnew


        

    
    def put_Sigma(self,Sigmaimp):
        """puts the impurity self energies for inequivalent atoms into the class.
        Finally, it should do the job of putting these inequiv. Sigmas at the correct points.
        """

        assert isinstance(Sigmaimp,list), "Sigmaimp has to be a list of Sigmas for the correlated shells, even if it is of length 1!"
        assert len(Sigmaimp)==self.N_inequiv_corr_shells, "give exactly one Sigma for each inequivalent corr. shell!"

       
        # init self.Sigmaimp:
        if (Sigmaimp[0].Note=='ReFreq'):
            # Real frequency Sigma:
            self.Sigmaimp = [ GF( Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = Sigmaimp[0].mesh)) for a,al in self.GFStruct_corr[i] ],
                                  Copy = False) for i in xrange(self.N_corr_shells) ]
        else:
            # Imaginary frequency Sigma:
            self.Sigmaimp = [ GF( Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = Sigmaimp[0].mesh)) for a,al in self.GFStruct_corr[i] ],
                                  Copy = False) for i in xrange(self.N_corr_shells) ]
                
        # transform the CTQMC blocks to the full matrix:
        for icrsh in xrange(self.N_corr_shells):
            s = self.shellmap[icrsh]    # s is the index of the inequivalent shell corresponding to icrsh
            for ibl in range(len(self.GFStruct_Solver[s])):
                for i in range(len(self.GFStruct_Solver[s][ibl][1])):
                    for j in range(len(self.GFStruct_Solver[s][ibl][1])):
                        bl   = self.GFStruct_Solver[s][ibl][0]
                        ind1 = self.GFStruct_Solver[s][ibl][1][i]
                        ind2 = self.GFStruct_Solver[s][ibl][1][j]
                        self.Sigmaimp[icrsh][self.mapinv[s][bl]][ind1,ind2] <<= Sigmaimp[s][bl][ind1,ind2]
                        
        #---------------------------------------------
        # double counting:
        # if the matrix dm is from the impurity problem(local coord. system) it has to be done here:
        #if hasattr(self,"dc_imp"): 
            #print "using double counting correction:"
            #print self.dc_imp
        #    for icrsh in xrange(self.N_corr_shells):
        #        for bl,gf in self.Sigmaimp[icrsh]: self.Sigmaimp[icrsh][bl] -= self.dc_imp[icrsh][bl]
        #---------------------------------------------
            
        

        # rotation from local to global coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in self.Sigmaimp[icrsh]: self.Sigmaimp[icrsh][sig] <<= self.rotloc(icrsh,gf,direction='toGlobal')



    def add_DC(self):

        sres = [s.copy() for s in self.Sigmaimp]
        if hasattr(self,"dc_imp"): 
            for icrsh in xrange(self.N_corr_shells):
                for bl,gf in sres[icrsh]: sres[icrsh][bl] -= self.dc_imp[icrsh][bl]
        #else:
        #    MPI.report("WARNING: No DC term set!!!")
            
        return sres 

               

    def set_mu(self,mu):
        """Sets a new chemical potential"""
        self.Chemical_Potential = mu
        #print "Chemical potential in SumK set to ",mu



    def sorts_of_atoms(self,lst):
        """
        This routine should determine the number of sorts in the double list lst
        """
        sortlst = [ lst[i][1] for i in xrange(len(lst)) ]
        sortlst.sort()
        sorts = 1
        for i in xrange(len(sortlst)-1):
            if sortlst[i+1]>sortlst[i]: sorts += 1

        return sorts



    def number_of_atoms(self,lst):
        """
        This routine should determine the number of atoms in the double list lst
        """
        atomlst = [ lst[i][0] for i in xrange(len(lst)) ]
        atomlst.sort()
        atoms = 1
        for i in xrange(len(atomlst)-1):
            if atomlst[i+1]>atomlst[i]: atoms += 1

        return atoms



    def inequiv_shells(self,lst):
        """
        The number of inequivalent shells is calculated from lst, and a mapping is given as
        map(i_corr_shells) = i_inequiv_corr_shells
        invmap(i_inequiv_corr_shells) = i_corr_shells
        in order to put the Self energies to all equivalent shells, and for extracting Gloc
        """

        tmp = []
        self.shellmap = [0 for i in range(len(lst))]
        self.invshellmap = [0]
        self.N_inequiv_corr_shells = 1
        tmp.append( lst[0][1:3] )
        
        if (len(lst)>1):
            for i in range(len(lst)-1):
               
                fnd = False
                for j in range(self.N_inequiv_corr_shells):
                    if (tmp[j]==lst[i+1][1:3]):
                        fnd = True
                        self.shellmap[i+1] = j
                if (fnd==False):
                    self.shellmap[i+1] = self.N_inequiv_corr_shells
                    self.N_inequiv_corr_shells += 1
                    tmp.append( lst[i+1][1:3] )
                    self.invshellmap.append(i+1)
                                

      
    def total_density(self, mu):
        """
        Calculates the total charge for the energy window for a given mu. Since in general N_Orbitals depends on k, 
        the calculation is done in the following order:
        G_aa'(k,iw) -> n(k) = Tr G_aa'(k,iw) -> sum_k n_k 
        
        mu: chemical potential
        
        If Projections do not depend on k, store the sum over k in self.__G in order to use it for the extraction of Gloc
        The calculation is done in the global coordinate system, if distinction is made between local/global!
        """

        
        if self.k_dep_projection == 0:
            assert 0,"In total_density, k-independent calculation has to be rewritten (spin-indices)!!"

            BS = range(self.N_Orbitals[0])
            self.__G = GF(Name_Block_Generator = [ (s,GFBloc_ImFreq(Indices = BS, Mesh = self.Sigmaimp[0].mesh)) for s in ['up','down'] ],Copy = False)
            self.__G.zero()
            tmp2 = self.__G.copy()       # initialise tmp2 for the k-sum already here
                                   
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in tmp2: 
                    tmp2[sig].from_L_G_R(self.Proj_Mat[icrsh].conjugate().transpose(),self.Sigmaimp[icrsh][sig],self.Proj_Mat[icrsh])  # upfolding this part
                self.__G += tmp2      # adding to the upfolded Sigma, use __G as temporary storage
                
            
            # now set things before starting the loop:
            mupat = numpy.identity(self.N_Orbitals[0],numpy.complex_)
            mupat *= mu

            tmp = tmp2.copy()             # initialise tmp with the correct size
            tmp <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
           
            tmp -= self.__G               # substract sigma, it does not depend on k!
           
            dens = 0.0
            self.__G.zero()               # now use __G to store the sum over k
           
            for ik in xrange(self.Nk):
                tmp2 <<= tmp
                tmp2 -= tmp2.NBlocks * [ self.Hopping[ik] - mupat ]       # this is for spin independent H(k) !!
                tmp2.invert()
                tmp2 *= self.BZ_weights[ik]

                self.__G += tmp2
                dens += tmp2.total_density()

            self.__G_calculated = True                                      # set a flag so that we know we have already done the sum over k
           
        else:
            
            # here we do not store __G, since for Gloc the projection is INSIDE the sum over k, which makes evaluation slower. 
            # The calculation of the density is done far more often than the calculation of Gloc!

            ntoi = self.names_to_ind[self.SO]
            bln = self.blocnames[self.SO]

            dens = 0.0

            BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
            GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
            S = GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct ],Copy = False)
            mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln]   # construct mupat
            for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu

            ikarray=numpy.array(range(self.Nk))
            stmp = self.add_DC()

            for ik in MPI.slice_array(ikarray):
                GFsize = [ gf.N1 for sig,gf in S]
                unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                       for ib in range(self.NspinblocsGF[self.SO]) ] )
  
                #if (self.N_Orbitals[ik]!=GFsize):
                if (not unchangedsize):
                    BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                    # construct the upfolded Sigma, if #bands changed
                    GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                    S = GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct ],Copy = False)
                    mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                    for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu
                                   

                S <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
                M = copy.deepcopy(mupat)
                for ibl in range(self.NspinblocsGF[self.SO]): 
                    ind = ntoi[bln[ibl]]
                    M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
                S -= M   

                tmp = S.copy()    # init temporary storage
                for icrsh in xrange(self.N_corr_shells):
                    for sig,gf in tmp: tmp[sig] <<= self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)
                    S -= tmp      # adding to the upfolded Sigma
                    
                S.invert()

                dens += self.BZ_weights[ik] * S.total_density()
                
            # collect data from MPI:
            dens = MPI.all_reduce(MPI.world,dens,lambda x,y : x+y)
            MPI.barrier()
                
        return dens


    def find_mu(self, precision = 0.01):
        """Searches for mu in order to give the desired charge
        A desired precision can be specified in precision."""

        F = lambda mu : self.total_density(mu = mu)

        Dens_rel = self.Density_Required - self.charge_below

        
        self.Chemical_Potential = Dichotomy.Dichotomy(Function = F,
                                         xinit = self.Chemical_Potential, yvalue = Dens_rel,
                                         Precision_on_y = precision, Delta_x=0.5,
                                         MaxNbreLoop = 100, xname="Chemical_Potential", yname= "Total Density",
                                         verbosity = 3)[0]

        return self.Chemical_Potential



    def find_mu_nonint(self, densreq, orb = None, Beta = 40, precision = 0.01):

        def F(mu):
            #gnonint = self.nonint_G(Beta=Beta,mu=mu)
            gnonint = self.extract_Gloc(mu=mu,withSigma=False)

            if (orb is None):
                dens = 0.0
                for ish in range(self.N_inequiv_corr_shells):
                    dens += gnonint[ish].total_density()    
            else:
                dens = gnonint[orb].total_density()
                
            return dens
            
        #F = lambda mu : self.nonint_G(Beta=Beta,mu=mu)[orb].total_density()

        self.Chemical_Potential = Dichotomy.Dichotomy(Function = F,
                                      xinit = self.Chemical_Potential, yvalue = densreq,
                                      Precision_on_y = precision, Delta_x=0.5,
                                      MaxNbreLoop = 100, xname="Chemical_Potential", yname= "Local Density",
                                      verbosity = 3)[0]

        return self.Chemical_Potential



    def extract_Gloc(self, mu=None, withSigma = True):
        """ 
        extracts the local downfolded Green function at the chemical potential of the class.
        At the end, the local G is rotated from the gloabl coordinate system to the local system.
        if withSigma = False: Sigma is not included => non-interacting local GF
        """

        if (mu is None): mu = self.Chemical_Potential

        
        if self.k_dep_projection == 0:

            if (self.__G_calculated==False):
                # the sum over k has not been done, density not calculated
                dens = self.total_density(mu=mu)

            Gloc = [ self.Sigmaimp[icrsh].copy() for icrsh in xrange(self.N_corr_shells) ]   # Gloc is the list of Gloc^R, be sure that it is a copy construction!
                                                                                               # this list will be returned  
            for icrsh in xrange(self.N_corr_shells): Gloc[icrsh].zero()                        # initialize to zero

            tmp = Gloc[0].copy()                                               # init temporay storage
            for icrsh in xrange(self.N_corr_shells):                           # downfolding G, out of k loop!
                for sig,gf in self.__G: 
                    tmp[sig].from_L_G_R(self.Proj_Mat[icrsh],self.__G[sig],self.Proj_Mat[icrsh].conjugate().transpose())
                Gloc[icrsh] += tmp          # += makes a copy and adds it to Gloc[icrsh]
                        
                                    
        else:

            ntoi = self.names_to_ind[self.SO]
            bln = self.blocnames[self.SO]

            # initialisation:
            BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
            GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
            S = GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct ],Copy = False)
            mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln]   # construct mupat
            for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu
            
            Gloc = [ self.Sigmaimp[icrsh].copy() for icrsh in xrange(self.N_corr_shells) ]   # this list will be returned  
            for icrsh in xrange(self.N_corr_shells): Gloc[icrsh].zero()                # initialize to zero

            ikarray=numpy.array(range(self.Nk))
            stmp = self.add_DC()

            for ik in MPI.slice_array(ikarray):
                GFsize = [ gf.N1 for sig,gf in S]  
                unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                       for ib in range(self.NspinblocsGF[self.SO]) ] )

                #if (self.N_Orbitals[ik]!=GFsize):
                if (not unchangedsize):
                    BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                    # construct the upfolded Sigma, if #bands changed
                    GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                    S = GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct ],Copy = False)
                    mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                    for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu
              
                S <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
                M = copy.deepcopy(mupat)
                for ibl in range(self.NspinblocsGF[self.SO]): 
                    ind = ntoi[bln[ibl]]
                    M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
                S -= M    

                # for tests without Sigma, skip following part:
                if (withSigma):
                    tmp = S.copy()    # init temporary storage
                    for icrsh in xrange(self.N_corr_shells):
                        for sig,gf in tmp: tmp[sig] <<= self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)
                        S -= tmp      # adding to the upfolded Sigma
             
                S.invert()
                S *= self.BZ_weights[ik]

                
                for icrsh in xrange(self.N_corr_shells):
                    tmp = Gloc[icrsh].copy()                  # init temporary storage
                    for sig,gf in tmp: tmp[sig] <<= self.downfold(ik,icrsh,sig,S[sig],gf)
                    Gloc[icrsh] += tmp

            #collect data from MPI:
            for icrsh in xrange(self.N_corr_shells):
                Gloc[icrsh] <<= MPI.all_reduce(MPI.world,Gloc[icrsh],lambda x,y : x+y)
            MPI.barrier()

  
        # Gloc[:] is now the sum over k projected to the local orbitals.
        # here comes the symmetrisation, if needed:   
        if (self.symm_op!=0): Gloc = self.Symm_corr.symmetrise(Gloc)
        
        # Gloc is rotated to the local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in Gloc[icrsh]: Gloc[icrsh][sig] <<= self.rotloc(icrsh,gf,direction='toLocal')

        # transform to CTQMC blocks:
        Glocret = [ GF( Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = Gloc[0].mesh)) for a,al in self.GFStruct_Solver[i] ],
                        Copy = False) for i in xrange(self.N_inequiv_corr_shells)  ]
        for ish in xrange(self.N_inequiv_corr_shells):
            for ibl in range(len(self.GFStruct_Solver[ish])):
                for i in range(len(self.GFStruct_Solver[ish][ibl][1])):
                    for j in range(len(self.GFStruct_Solver[ish][ibl][1])):
                        bl   = self.GFStruct_Solver[ish][ibl][0]
                        ind1 = self.GFStruct_Solver[ish][ibl][1][i]
                        ind2 = self.GFStruct_Solver[ish][ibl][1][j]
                        Glocret[ish][bl][ind1,ind2] <<= Gloc[self.invshellmap[ish]][self.mapinv[ish][bl]][ind1,ind2]


        # return only the inequivalent shells:
        return Glocret
   

    def calc_DensityCorrection(self,Filename = 'densmat.dat'):
        """ Calculates the density correction in order to feed it back to the DFT calculations."""

        
        assert (type(Filename)==StringType), "Filename has to be a string!"

        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        # Set up deltaN:
        deltaN = {}
        for ib in bln:
            deltaN[ib] = [ numpy.zeros( [self.N_Orbitals[ik][ntoi[ib]],self.N_Orbitals[ik][ntoi[ib]]], numpy.complex_) for ik in range(self.Nk)]

        # initialisation:
        BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
        GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
        S = GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct ],Copy = False)
        mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln]   # construct mupat
        for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential
        
               
        ikarray=numpy.array(range(self.Nk))
        #ikarray=numpy.array([0])
        dens = 0.0
        stmp = self.add_DC()

        for ik in MPI.slice_array(ikarray):
            GFsize = [ gf.N1 for sig,gf in S]  
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                   for ib in range(self.NspinblocsGF[self.SO]) ] )

            if (not unchangedsize):
                BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                # construct the upfolded Sigma, if #bands changed
                GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                S = GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct ],Copy = False)
                mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential
              
            S <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
            M = copy.deepcopy(mupat)
            for ibl in range(self.NspinblocsGF[self.SO]): 
                ind = ntoi[bln[ibl]]
                M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
            S -= M    
            # S is now the non-interacting Green function
            
            Gnonint = S.copy()
            Gnonint.invert()

            # for tests without Sigma, skip following part:
            tmp = S.copy()    # init temporary storage
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in tmp: tmp[sig] <<= self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)
                S -= tmp      # adding to the upfolded Sigma
             
            S.invert()

            #Gnonint.save('Gnonint')
            #S.save('Gint')
            
            # first possibility:
            # calculate (G_int - G_nonint):
            #S -= Gnonint

            # second possibility:
            # Gnonint*Sigma*Gint:
            #S <<= (Gnonint*tmp)*S

            # third possibility:
            #for sig,f in S:
            #    deltaN[sig][ik] = S[sig].density() - Gnonint[sig].density()

            for sig,g in S:
                deltaN[sig][ik] = S[sig].density()
            dens += self.BZ_weights[ik] * S.total_density()
            
                

        #put MPI Barrier:
        for sig in deltaN:
            for ik in range(self.Nk):
                deltaN[sig][ik] = MPI.all_reduce(MPI.world,deltaN[sig][ik],lambda x,y : x+y)
        dens = MPI.all_reduce(MPI.world,dens,lambda x,y : x+y)
        MPI.barrier()

       
        # now save to file:
        if (MPI.IS_MASTER_NODE()):
            if (self.SP==0):
                f=open(Filename,'w')
            else:
                f=open(Filename+'up','w')
                f1=open(Filename+'dn','w')
            # write chemical potential (in Rydberg):
            f.write("%.14f\n"%(self.Chemical_Potential/self.EnergyUnit))
            if (self.SP!=0): f1.write("%.14f\n"%(self.Chemical_Potential/self.EnergyUnit))
            # write beta in ryderg-1
            f.write("%.14f\n"%(S.Beta*self.EnergyUnit))
            if (self.SP!=0): f1.write("%.14f\n"%(S.Beta*self.EnergyUnit))
            if (self.SP==0):
                for ik in range(self.Nk):
                    f.write("%s\n"%self.N_Orbitals[ik][0])
                    for inu in range(self.N_Orbitals[ik][0]):
                        for imu in range(self.N_Orbitals[ik][0]):
                            f.write("%.14f  %.14f "%(deltaN['up'][ik][inu,imu].real,deltaN['up'][ik][inu,imu].imag))
                        f.write("\n")
                    f.write("\n")
                f.close()
            else:
                 for ik in range(self.Nk):
                     f.write("%s\n"%self.N_Orbitals[ik][0])
                     for inu in range(self.N_Orbitals[ik][0]):
                         for imu in range(self.N_Orbitals[ik][0]):
                             f.write("%.14f  %.14f "%(deltaN['up'][ik][inu,imu].real,deltaN['up'][ik][inu,imu].imag))
                         f.write("\n")
                     f.write("\n")
                 f.close()
                 for ik in range(self.Nk):
                     f1.write("%s\n"%self.N_Orbitals[ik][1])
                     for inu in range(self.N_Orbitals[ik][1]):
                         for imu in range(self.N_Orbitals[ik][1]):
                             f1.write("%.14f  %.14f "%(deltaN['down'][ik][inu,imu].real,deltaN['down'][ik][inu,imu].imag))
                         f1.write("\n")
                     f1.write("\n")
                 f1.close()


        return deltaN, dens
