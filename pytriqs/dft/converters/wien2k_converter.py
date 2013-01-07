
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
import numpy
from pytriqs.base.archive import *
import pytriqs.base.utility.mpi as mpi
import string


def read_fortran_file (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
	    yield string.atof(x)



class Wien2kConverter:
    """
    Conversion from Wien2k output to an hdf5 file, that can be used as input for the SumkLDA class.
    """

    def __init__(self, filename, lda_subgrp = 'SumK_LDA', summ_subgrp = 'SymmCorr', repacking = False):
        """
        Init of the class. Variable filename gives the root of all filenames, e.g. case.ctqmcout, case.h5, and so
        on. 
        """

        assert type(filename)==StringType,"LDA_file must be a filename"
        self.HDFfile = filename+'.h5'
        self.LDA_file = filename+'.ctqmcout'
        self.Symm_file = filename+'.symqmc'
        self.Parproj_file = filename+'.parproj'
        self.Symmpar_file = filename+'.sympar'
        self.Band_file = filename+'.outband'
        self.lda_subgrp = lda_subgrp
        self.summ_subgrp = summ_subgrp

        # Checks if h5 file is there and repacks it if wanted:
        import os.path
        if (os.path.exists(self.HDFfile) and repacking):
            self.__repack()
        
        

    def convert_dmft_input(self):
        """
        Reads the input files, and stores the data in the HDFfile
        """
        
                   
        if not (mpi.is_master_node()): return # do it only on master:
        mpi.report("Reading input from %s..."%self.LDA_file)

        # Read and write only on Master!!!
        # R is a generator : each R.Next() will return the next number in the file
        R = read_fortran_file(self.LDA_file)
        try:
            EnergyUnit = R.next()                         # read the energy convertion factor
            Nk = int(R.next())                            # read the number of k points
            k_dep_projection = 1                          
            SP = int(R.next())                            # flag for spin-polarised calculation
            SO = int(R.next())                            # flag for spin-orbit calculation
            charge_below = R.next()                       # total charge below energy window
            Density_Required = R.next()                   # total density required, for setting the chemical potential
            symm_op = 1                                   # Use symmetry groups for the k-sum

            # the information on the non-correlated shells is not important here, maybe skip:
            N_shells = int(R.next())                      # number of shells (e.g. Fe d, As p, O p) in the unit cell, 
                                                               # corresponds to index R in formulas
            # now read the information about the shells:
            shells = [ [ int(R.next()) for i in range(4) ] for icrsh in range(N_shells) ]    # reads iatom, sort, l, dim

            N_corr_shells = int(R.next())                 # number of corr. shells (e.g. Fe d, Ce f) in the unit cell, 
                                                          # corresponds to index R in formulas
            # now read the information about the shells:
            corr_shells = [ [ int(R.next()) for i in range(6) ] for icrsh in range(N_corr_shells) ]    # reads iatom, sort, l, dim, SO flag, irep

            self.inequiv_shells(corr_shells)              # determine the number of inequivalent correlated shells, has to be known for further reading...


            use_rotations = 1
            rotmat = [numpy.identity(corr_shells[icrsh][3],numpy.complex_) for icrsh in xrange(N_corr_shells)]
           
            # read the matrices
            rotmat_timeinv = [0 for i in range(N_corr_shells)]

            for icrsh in xrange(N_corr_shells):
                for i in xrange(corr_shells[icrsh][3]):    # read real part:
                    for j in xrange(corr_shells[icrsh][3]):
                        rotmat[icrsh][i,j] = R.next()
                for i in xrange(corr_shells[icrsh][3]):    # read imaginary part:
                    for j in xrange(corr_shells[icrsh][3]):
                        rotmat[icrsh][i,j] += 1j * R.next()

                if (SP==1):             # read time inversion flag:
                    rotmat_timeinv[icrsh] = int(R.next())
                    
                  
            
            # Read here the infos for the transformation of the basis:
            Nreps = [1 for i in range(self.N_inequiv_corr_shells)]
            dim_reps = [0 for i in range(self.N_inequiv_corr_shells)]
            T = []
            for icrsh in range(self.N_inequiv_corr_shells):
                Nreps[icrsh] = int(R.next())   # number of representatives ("subsets"), e.g. t2g and eg
                dim_reps[icrsh] = [int(R.next()) for i in range(Nreps[icrsh])]   # dimensions of the subsets
            
            # The transformation matrix:
            # it is of dimension 2l+1, if no SO, and 2*(2l+1) with SO!!
            #T = []
            #for ish in xrange(self.N_inequiv_corr_shells):
                ll = 2*corr_shells[self.invshellmap[icrsh]][2]+1
                lmax = ll * (corr_shells[self.invshellmap[icrsh]][4] + 1)
                T.append(numpy.zeros([lmax,lmax],numpy.complex_))
                
                # now read it from file:
                for i in xrange(lmax):
                    for j in xrange(lmax):
                        T[icrsh][i,j] = R.next()
                for i in xrange(lmax):
                    for j in xrange(lmax):
                        T[icrsh][i,j] += 1j * R.next()

    
            # Spin blocks to be read:
            Nspinblocs = SP + 1 - SO   # number of spins to read for Norbs and Ham, NOT Projectors
                 
        
            # read the list of N_Orbitals for all k points
            N_Orbitals = [ [0 for isp in range(Nspinblocs)] for ik in xrange(Nk)]
            for isp in range(Nspinblocs):
                for ik in xrange(Nk):
                    N_Orbitals[ik][isp] = int(R.next())
            #print N_Orbitals

            # Initialise the projectors:
            Proj_Mat = [ [ [numpy.zeros([corr_shells[icrsh][3], N_Orbitals[ik][isp]], numpy.complex_) 
                            for icrsh in range (N_corr_shells)] 
                           for isp in range(Nspinblocs)] 
                         for ik in range(Nk) ]

            # Read the projectors from the file:
            for ik in xrange(Nk):
                for icrsh in range(N_corr_shells):
                    no = corr_shells[icrsh][3]
                    # first Real part for BOTH spins, due to conventions in dmftproj:
                    for isp in range(Nspinblocs):
                        for i in xrange(no):
                            for j in xrange(N_Orbitals[ik][isp]):
                                Proj_Mat[ik][isp][icrsh][i,j] = R.next()
                    # now Imag part:
                    for isp in range(Nspinblocs):
                        for i in xrange(no):
                            for j in xrange(N_Orbitals[ik][isp]):
                                Proj_Mat[ik][isp][icrsh][i,j] += 1j * R.next()
            
          
            # now define the arrays for weights and hopping ...
            BZ_weights = numpy.ones([Nk],numpy.float_)/ float(Nk)  # w(k_index),  default normalisation 
            Hopping = [ [numpy.zeros([N_Orbitals[ik][isp],N_Orbitals[ik][isp]],numpy.complex_) 
                         for isp in range(Nspinblocs)] for ik in xrange(Nk) ]

                            
            # weights in the file
            for ik in xrange(Nk) : BZ_weights[ik] = R.next()         
                
            # if the sum over spins is in the weights, take it out again!!
            sm = sum(BZ_weights)
            BZ_weights[:] /= sm 
	    
            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian!!!!
            for isp in range(Nspinblocs):
                for ik in xrange(Nk) :
                    no = N_Orbitals[ik][isp]
                    for i in xrange(no):
                        Hopping[ik][isp][i,i] = R.next() * EnergyUnit
            
            #keep some things that we need for reading parproj:
            self.N_shells = N_shells
            self.shells = shells
            self.N_corr_shells = N_corr_shells
            self.corr_shells = corr_shells
            self.Nspinblocs = Nspinblocs
            self.N_Orbitals = N_Orbitals
            self.Nk = Nk
            self.SO = SO
            self.SP = SP
            self.EnergyUnit = EnergyUnit
        except StopIteration : # a more explicit error if the file is corrupted.
            raise "SumkLDA : reading file HMLT_file failed!"

        R.close()
        
        #print Proj_Mat[0]

        #-----------------------------------------
        # Store the input into HDF5:
        ar = HDFArchive(self.HDFfile,'a')
        if not (self.lda_subgrp in ar): ar.create_group(self.lda_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created.
        # If it exists, the data is overwritten!!!
        
        ar[self.lda_subgrp]['EnergyUnit'] = EnergyUnit
        ar[self.lda_subgrp]['Nk'] = Nk
        ar[self.lda_subgrp]['k_dep_projection'] = k_dep_projection
        ar[self.lda_subgrp]['SP'] = SP
        ar[self.lda_subgrp]['SO'] = SO
        ar[self.lda_subgrp]['charge_below'] = charge_below
        ar[self.lda_subgrp]['Density_Required'] = Density_Required
        ar[self.lda_subgrp]['symm_op'] = symm_op
        ar[self.lda_subgrp]['N_shells'] = N_shells
        ar[self.lda_subgrp]['shells'] = shells
        ar[self.lda_subgrp]['N_corr_shells'] = N_corr_shells
        ar[self.lda_subgrp]['corr_shells'] = corr_shells
        ar[self.lda_subgrp]['use_rotations'] = use_rotations
        ar[self.lda_subgrp]['rotmat'] = rotmat
        ar[self.lda_subgrp]['rotmat_timeinv'] = rotmat_timeinv
        ar[self.lda_subgrp]['Nreps'] = Nreps
        ar[self.lda_subgrp]['dim_reps'] = dim_reps
        ar[self.lda_subgrp]['T'] = T
        ar[self.lda_subgrp]['N_Orbitals'] = N_Orbitals
        ar[self.lda_subgrp]['Proj_Mat'] = Proj_Mat
        ar[self.lda_subgrp]['BZ_weights'] = BZ_weights
        ar[self.lda_subgrp]['Hopping'] = Hopping
        
        del ar
              
        # Symmetries are used, 
        # Now do the symmetries for correlated orbitals:
        self.read_symmetry_input(orbits=corr_shells,symm_file=self.Symm_file,summ_subgrp=self.summ_subgrp,SO=SO,SP=SP)


    def convert_parproj_input(self, par_proj_subgrp='SumK_LDA_ParProj', symm_par_subgrp='SymmPar'):
        """
        Reads the input for the partial charges projectors from case.parproj, and stores it in the symm_par_subgrp
        group in the HDF5.
        """

        if not (mpi.is_master_node()): return

        self.par_proj_subgrp = par_proj_subgrp
        self.symm_par_subgrp = symm_par_subgrp

        mpi.report("Reading parproj input from %s..."%self.Parproj_file)

        Dens_Mat_below = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.N_shells)] 
                           for isp in range(self.Nspinblocs) ]

        R = read_fortran_file(self.Parproj_file)
        #try:

        N_parproj = [int(R.next()) for i in range(self.N_shells)]
                
        # Initialise P, here a double list of matrices:
        Proj_Mat_pc = [ [ [ [numpy.zeros([self.shells[ish][3], self.N_Orbitals[ik][isp]], numpy.complex_) 
                             for ir in range(N_parproj[ish])]
                            for ish in range (self.N_shells) ]
                          for isp in range(self.Nspinblocs) ]
                        for ik in range(self.Nk) ]

        rotmat_all = [numpy.identity(self.shells[ish][3],numpy.complex_) for ish in xrange(self.N_shells)]
        rotmat_all_timeinv = [0 for i in range(self.N_shells)]

        for ish in range(self.N_shells):
            #print ish   
            # read first the projectors for this orbital:
            for ik in xrange(self.Nk):
                for ir in range(N_parproj[ish]):
                    for isp in range(self.Nspinblocs):
                                    
                        for i in xrange(self.shells[ish][3]):    # read real part:
                            for j in xrange(self.N_Orbitals[ik][isp]):
                                Proj_Mat_pc[ik][isp][ish][ir][i,j] = R.next()
                            
                    for isp in range(self.Nspinblocs):
                        for i in xrange(self.shells[ish][3]):    # read imaginary part:
                            for j in xrange(self.N_Orbitals[ik][isp]):
                                Proj_Mat_pc[ik][isp][ish][ir][i,j] += 1j * R.next()
                                        
                    
            # now read the Density Matrix for this orbital below the energy window:
            for isp in range(self.Nspinblocs):
                for i in xrange(self.shells[ish][3]):    # read real part:
                    for j in xrange(self.shells[ish][3]):
                        Dens_Mat_below[isp][ish][i,j] = R.next()
            for isp in range(self.Nspinblocs):
                for i in xrange(self.shells[ish][3]):    # read imaginary part:
                    for j in xrange(self.shells[ish][3]):
                        Dens_Mat_below[isp][ish][i,j] += 1j * R.next()
                if (self.SP==0): Dens_Mat_below[isp][ish] /= 2.0

            # Global -> local rotation matrix for this shell:
            for i in xrange(self.shells[ish][3]):    # read real part:
                for j in xrange(self.shells[ish][3]):
                    rotmat_all[ish][i,j] = R.next()
            for i in xrange(self.shells[ish][3]):    # read imaginary part:
                for j in xrange(self.shells[ish][3]):
                    rotmat_all[ish][i,j] += 1j * R.next()
                    
            #print Dens_Mat_below[0][ish],Dens_Mat_below[1][ish]
            
            if (self.SP):
                rotmat_all_timeinv[ish] = int(R.next())

        #except StopIteration : # a more explicit error if the file is corrupted.
        #    raise "Wien2kConverter: reading file for Projectors failed!"
        R.close()

        #-----------------------------------------
        # Store the input into HDF5:
        ar = HDFArchive(self.HDFfile,'a')
        if not (self.par_proj_subgrp in ar): ar.create_group(self.par_proj_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created.
        # If it exists, the data is overwritten!!!
        thingstowrite = ['Dens_Mat_below','N_parproj','Proj_Mat_pc','rotmat_all','rotmat_all_timeinv']
        for it in thingstowrite: exec "ar['%s']['%s'] = %s"%(self.par_proj_subgrp,it,it)
        del ar

        # Symmetries are used, 
        # Now do the symmetries for all orbitals:
        self.read_symmetry_input(orbits=self.shells,symm_file=self.Symmpar_file,summ_subgrp=self.symm_par_subgrp,SO=self.SO,SP=self.SP)


    def convert_bands_input(self, bands_subgrp = 'SumK_LDA_Bands'):
        """
        Converts the input for momentum resolved spectral functions, and stores it in bands_subgrp in the
        HDF5.
        """

        if not (mpi.is_master_node()): return

        self.bands_subgrp = bands_subgrp
        mpi.report("Reading bands input from %s..."%self.Band_file)

        R = read_fortran_file(self.Band_file)
        try:
            Nk = int(R.next())

            # read the list of N_Orbitals for all k points
            N_Orbitals = [ [0 for isp in range(self.Nspinblocs)] for ik in xrange(Nk)]
            for isp in range(self.Nspinblocs):
                for ik in xrange(Nk):
                    N_Orbitals[ik][isp] = int(R.next())

            # Initialise the projectors:
            Proj_Mat = [ [ [numpy.zeros([self.corr_shells[icrsh][3], N_Orbitals[ik][isp]], numpy.complex_) 
                            for icrsh in range (self.N_corr_shells)] 
                           for isp in range(self.Nspinblocs)] 
                         for ik in range(Nk) ]

            # Read the projectors from the file:
            for ik in xrange(Nk):
                for icrsh in range(self.N_corr_shells):
                    no = self.corr_shells[icrsh][3]
                    # first Real part for BOTH spins, due to conventions in dmftproj:
                    for isp in range(self.Nspinblocs):
                        for i in xrange(no):
                            for j in xrange(N_Orbitals[ik][isp]):
                                Proj_Mat[ik][isp][icrsh][i,j] = R.next()
                    # now Imag part:
                    for isp in range(self.Nspinblocs):
                        for i in xrange(no):
                            for j in xrange(N_Orbitals[ik][isp]):
                                Proj_Mat[ik][isp][icrsh][i,j] += 1j * R.next()

            Hopping = [ [numpy.zeros([N_Orbitals[ik][isp],N_Orbitals[ik][isp]],numpy.complex_) 
                         for isp in range(self.Nspinblocs)] for ik in xrange(Nk) ]
         	    
            # Grab the H
            # we use now the convention of a DIAGONAL Hamiltonian!!!!
            for isp in range(self.Nspinblocs):
                for ik in xrange(Nk) :
                    no = N_Orbitals[ik][isp]
                    for i in xrange(no):
                        Hopping[ik][isp][i,i] = R.next() * self.EnergyUnit

            # now read the partial projectors:
            N_parproj = [int(R.next()) for i in range(self.N_shells)]
            # Initialise P, here a double list of matrices:
            Proj_Mat_pc = [ [ [ [numpy.zeros([self.shells[ish][3], N_Orbitals[ik][isp]], numpy.complex_) 
                                 for ir in range(N_parproj[ish])]
                                for ish in range (self.N_shells) ]
                              for isp in range(self.Nspinblocs) ]
                            for ik in range(Nk) ]


            for ish in range(self.N_shells):
               
                for ik in xrange(Nk):
                    for ir in range(N_parproj[ish]):
                        for isp in range(self.Nspinblocs):
                                    
                            for i in xrange(self.shells[ish][3]):    # read real part:
                                for j in xrange(N_Orbitals[ik][isp]):
                                    Proj_Mat_pc[ik][isp][ish][ir][i,j] = R.next()
                            
                            for i in xrange(self.shells[ish][3]):    # read imaginary part:
                                for j in xrange(N_Orbitals[ik][isp]):
                                    Proj_Mat_pc[ik][isp][ish][ir][i,j] += 1j * R.next()

        except StopIteration : # a more explicit error if the file is corrupted.
            raise "SumkLDA : reading file HMLT_file failed!"

        R.close()
        # reading done!

        #-----------------------------------------
        # Store the input into HDF5:
        ar = HDFArchive(self.HDFfile,'a')
        if not (self.bands_subgrp in ar): ar.create_group(self.bands_subgrp) 
        # The subgroup containing the data. If it does not exist, it is created.
        # If it exists, the data is overwritten!!!
        thingstowrite = ['Nk','N_Orbitals','Proj_Mat','Hopping','N_parproj','Proj_Mat_pc']
        for it in thingstowrite: exec "ar['%s']['%s'] = %s"%(self.bands_subgrp,it,it)

        #ar[self.bands_subgrp]['Nk'] = Nk
        #ar[self.bands_subgrp]['N_Orbitals'] = N_Orbitals
        #ar[self.bands_subgrp]['Proj_Mat'] = Proj_Mat
        #self.Proj_Mat = Proj_Mat
        #self.N_Orbitals = N_Orbitals
        #self.Nk = Nk
        #self.Hopping = Hopping
        del ar
   




    def read_symmetry_input(self, orbits, symm_file, summ_subgrp, SO, SP):
        """
        Reads input for the symmetrisations from symm_file, which is case.sympar or case.symqmc.
        """

        if not (mpi.is_master_node()): return

        mpi.report("Reading symmetry input from %s..."%symm_file)

        N_orbits = len(orbits)
        R=read_fortran_file(symm_file)

        try:
            Ns = int(R.next())           # Number of symmetry operations
            Natoms = int(R.next())       # number of atoms involved
            perm = [ [int(R.next()) for i in xrange(Natoms)] for j in xrange(Ns) ]    # list of permutations of the atoms
            if SP: 
                timeinv = [ int(R.next()) for j in xrange(Ns) ]           # timeinversion for SO xoupling
            else:
                timeinv = [ 0 for j in xrange(Ns) ] 

            # Now read matrices:
            mat = []  
            for iNs in xrange(Ns):
                
                mat.append( [ numpy.zeros([orbits[orb][3], orbits[orb][3]],numpy.complex_) for orb in xrange(N_orbits) ] )
                for orb in range(N_orbits):
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat[iNs][orb][i,j] = R.next()            # real part
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat[iNs][orb][i,j] += 1j * R.next()      # imaginary part

            # determine the inequivalent shells:
            #SHOULD BE FINALLY REMOVED, PUT IT FOR ALL ORBITALS!!!!!
            #self.inequiv_shells(orbits)
            mat_tinv = [numpy.identity(orbits[orb][3],numpy.complex_)
                        for orb in range(N_orbits)]

            if ((SO==0) and (SP==0)):
                # here we need an additional time inversion operation, so read it:
                for orb in range(N_orbits):
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat_tinv[orb][i,j] = R.next()            # real part
                    for i in xrange(orbits[orb][3]):
                        for j in xrange(orbits[orb][3]):
                            mat_tinv[orb][i,j] += 1j * R.next()      # imaginary part
                


        except StopIteration : # a more explicit error if the file is corrupted.
	    raise "Symmetry : reading file failed!"
        
        R.close()

        # Save it to the HDF:
        ar=HDFArchive(self.HDFfile,'a')
        if not (summ_subgrp in ar): ar.create_group(summ_subgrp)
        thingstowrite = ['Ns','Natoms','perm','orbits','SO','SP','timeinv','mat','mat_tinv']
        for it in thingstowrite: exec "ar['%s']['%s'] = %s"%(summ_subgrp,it,it)
        del ar
        
        

    def __repack(self):
        """Calls the h5repack routine, in order to reduce the file size of the hdf5 archive.
           Should only be used BEFORE the first invokation of HDFArchive in the program, otherwise
           the hdf5 linking is broken!!!"""

        import subprocess

        if not (mpi.is_master_node()): return

        mpi.report("Repacking the file %s"%self.HDFfile)

        retcode = subprocess.call(["h5repack","-i%s"%self.HDFfile, "-otemphgfrt.h5"])
        if (retcode!=0):
            mpi.report("h5repack failed!")
        else:
            subprocess.call(["mv","-f","temphgfrt.h5","%s"%self.HDFfile])
            


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
                                
