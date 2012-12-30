
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

#from pytriqs.import myUtils 
from types import *
import numpy
import pytriqs.base.Utility.Dichotomy as Dichotomy
from pytriqs.base.GF_Local.GF import GF
from pytriqs.base.GF_Local.GFBloc_ImFreq import GFBloc_ImFreq
from pytriqs.base.GF_Local.GFBloc_ReFreq import GFBloc_ReFreq
from pytriqs.base.GF_Local.GFBloc_ImTime import GFBloc_ImTime
from pytriqs.base.GF_Local import GF_Initializers
from pytriqs.solvers.operators import *
from pytriqs.base.Utility.myUtils import Sum
import pytriqs.base.Utility.MPI as MPI

from pytriqs.Wien2k.Symmetry import *
from pytriqs.Wien2k.SumK_LDA_SO import SumK_LDA_SO

import string


def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
	    yield string.atof(x)

def Read_Fortran_File2(filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
            try:
                yield string.atof(x)
            except:
                yield x.replace('E','D')


class SumK_LDA_SOtools(SumK_LDA_SO):
    """Extends the SumK_LDA class with some tools for analysing the data."""




   
    

    def downfold_pc(self,ik,ir,ish,sig,gf_to_downfold,gf_inp):
        """Downfolding a block of the Greens function"""
        
        gf_downfolded = gf_inp.copy()
        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        gf_downfolded.from_L_G_R(self.Proj_Mat_pc[ik][isp][ish][ir],gf_to_downfold,self.Proj_Mat_pc[ik][isp][ish][ir].conjugate().transpose())  # downfolding G

        return gf_downfolded


    def rotloc_all(self,ish,gf_to_rotate,direction):
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

        assert ((direction=='toLocal')or(direction=='toGlobal')),"Give direction 'toLocal' or 'toGlobal' in rotloc!"

        gf_rotated = gf_to_rotate.copy()
        if (direction=='toLocal'):
            gf_rotated.from_L_G_R(self.rotmat_all[ish].transpose(),gf_to_rotate,self.rotmat_all[ish].conjugate())
        elif (direction=='toGlobal'):
            gf_rotated.from_L_G_R(self.rotmat_all[ish].conjugate(),gf_to_rotate,self.rotmat_all[ish].transpose())
    
        return gf_rotated



    def check_inputDOS(self, ommin, ommax, N_om, Beta=10, broadening=0.01):
        
       
        delta_om = (ommax-ommin)/(N_om-1)
            
        Mesh = numpy.zeros([N_om],numpy.float_)

        DOS = {}
        for bn in self.blocnames[self.SO]:
            DOS[bn] = numpy.zeros([N_om],numpy.float_)

        DOSproj     = [ {} for icrsh in range(self.N_inequiv_corr_shells) ]
        DOSproj_orb = [ {} for icrsh in range(self.N_inequiv_corr_shells) ]
        for icrsh in range(self.N_inequiv_corr_shells):
            for bn in self.blocnames[self.corr_shells[self.invshellmap[icrsh]][4]]:
                dl = self.corr_shells[self.invshellmap[icrsh]][3]*(1+self.corr_shells[self.invshellmap[icrsh]][4])
                DOSproj[icrsh][bn] = numpy.zeros([N_om],numpy.float_)
                DOSproj_orb[icrsh][bn] = numpy.zeros([dl,dl,N_om],numpy.float_)

        #for bn in self.blocnames[self.SO]
        #DOS = [ numpy.zeros([N_om],numpy.float_) for bn in self.blocnames[self.SO] ]
        #DOSproj = [ numpy.zeros([N_om],numpy.float_) for icrsh in range(self.N_inequiv_corr_shells) ]
        #DOSproj_orb = [ numpy.zeros([self.corr_shells[self.invshellmap[ish]][3],self.corr_shells[self.invshellmap[ish]][3],N_om ],numpy.float_) 
        #                for ish in range(self.N_inequiv_corr_shells) ]

        for i in range(N_om): Mesh[i] = ommin + delta_om * i

        # init:
        Gloc = []
        for icrsh in range(self.N_corr_shells):
            b_list = [a for a,al in self.GFStruct_corr[icrsh]]
            glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in self.GFStruct_corr[icrsh]]   
            Gloc.append(GF(NameList = b_list, BlockList = glist(),Copy=False))
        for icrsh in xrange(self.N_corr_shells): Gloc[icrsh].zero()                        # initialize to zero
            

        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
        GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
        a_list = [a for a,al in GFStruct] 
        glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in GFStruct]       # Indices for the upfolded G
        Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
        Gupf.zero()
        mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln] 
        for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential


        #GFStruct = [ (ib, BS[self.names_to_ind[ib]]) for ib in self.blocnames ]
        #a_list = [a for a,al in GFStruct]   
        #glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in GFStruct]       # Indices for the upfolded G
        #Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
        #Gupf.zero()
        #mupat = [numpy.identity(self.N_Orbitals[0][self.names_to_ind[bl]],numpy.complex_) for bl in self.blocnames] 
        #for ibl in range(len(self.blocnames)): mupat[ibl] *= self.Chemical_Potential

        #for ik in xrange(1):
        for ik in xrange(self.Nk):

            if (ik%20==0): print ik

            GFsize = [ gf.N1 for sig,gf in Gupf]  
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                   for ib in range(self.NspinblocsGF[self.SO]) ] )

            if (not unchangedsize):
                BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                a_list = [a for a,al in GFStruct]                                 
                glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in GFStruct]    
                Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
                Gupf.zero()
                mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= self.Chemical_Potential
            
                
            Gupf <<= GF_Initializers.A_Omega_Plus_B(A=1,B=1j*broadening)
            M = copy.deepcopy(mupat)
            for ibl in range(self.NspinblocsGF[self.SO]): 
                ind = ntoi[bln[ibl]]
                M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
          
            Gupf -= M
                
            Gupf.invert()
            Gupf *= self.BZ_weights[ik]
            #print ik, Gupf

            # non-projected DOS
            for iom in range(N_om): 
                for sig,gf in Gupf: 
                    asd = gf._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                    DOS[sig][iom] += asd
                
            for icrsh in xrange(self.N_corr_shells):
                tmp = Gloc[icrsh].copy()
                for sig,gf in tmp: tmp[sig] <<= self.downfold(ik,icrsh,sig,Gupf[sig],gf) # downfolding G
                Gloc[icrsh] += tmp

                            
        
        if (self.symm_op!=0): Gloc = self.Symm_corr.symmetrise(Gloc)

        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in Gloc[icrsh]: Gloc[icrsh][sig] <<= self.rotloc(icrsh,gf,direction='toLocal')
                
        # Gloc can now also be used to look at orbitally resolved quantities
        for ish in range(self.N_inequiv_corr_shells):
            for sig,gf in Gloc[self.invshellmap[ish]]: # loop over spins
                for iom in range(N_om): DOSproj[ish][sig][iom] += gf._data.array[:,:,iom].imag.trace()/(-3.1415926535) 

                DOSproj_orb[ish][sig][:,:,:] += gf._data.array[:,:,:].imag/(-3.1415926535)
       
        if ((self.SP==0) and (self.SO==0)):
            DOS['up'][:] *= 2.0
            for ish in range(self.N_inequiv_corr_shells):
                DOSproj[ish]['up'][:] *= 2.0
                DOSproj_orb[ish]['up'][:,:,:] *= 2.0
          

        # output:
        #if (self.Nspinblocs==1):
        if (MPI.IS_MASTER_NODE()):
            for ibn in range(self.Nspinblocs):
                bn = self.blocnames[self.SO][ibn]
                f=open('DOS%s.dat'%bn, 'w')
                for i in range(N_om): f.write("%s    %s\n"%(Mesh[i],DOS[bn][i]))
                f.close()  

                for ish in range(self.N_inequiv_corr_shells):
                    f=open('DOS%sproj%s.dat'%(bn,ish),'w')
                    for i in range(N_om): f.write("%s    %s\n"%(Mesh[i],DOSproj[ish][bn][i]))
                    f.close()  
            
                    for i in range(self.corr_shells[self.invshellmap[ish]][3]):
                        for j in range(i,self.corr_shells[self.invshellmap[ish]][3]):
                            Fname = 'DOS'+bn+'proj_'+str(ish)+'_'+str(i)+'_'+str(j)+'.dat'
                            f=open(Fname,'w')
                            for iom in range(N_om): f.write("%s    %s\n"%(Mesh[iom],DOSproj_orb[ish][bn][i,j,iom]))
                            f.close()
        


    def DOSpartial(self,broadening=0.01,Proj_filename='',Symm_file='Symmetry.dat'):
        """should calculate the orbital-resolved DOS"""

        assert hasattr(self,"Sigmaimp"), "Set Sigma First!!"

        if (self.SO==0):
            self.__readprojfile(Proj_filename=Proj_filename,Symm_file=Symm_file)
        else:
            MPI.report("WARNING: with SO we calculate only TOTAL DOS for the moment!")

        #GFStruct_proj has to be adapted for SO coupling!!
        #GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in self.blocnames ]  for i in xrange(self.N_shells) ]
	#Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct_proj[ish] ], Copy = False ) 
	#	 for ish in xrange(self.N_shells)]
	#for ish in range(self.N_shells): Gproj[ish].zero()

        mu = self.Chemical_Potential
        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        if (self.SO==0):
            GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in bln ]  for i in xrange(self.N_shells) ]
            Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct_proj[ish] ], Copy = False ) 
                     for ish in xrange(self.N_shells)]
            for ish in range(self.N_shells): Gproj[ish].zero()

                
        if self.k_dep_projection == 0:
            assert 0,"Not implemented!!"	
        else:
             BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
             GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
             a_list = [a for a,al in GFStruct]   
             glist = lambda : [ GFBloc_ReFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]  
             S = GF(NameList = a_list, BlockList = glist(),Copy=False)
             S.zero()
             mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln] 
             for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu

             Msh = [x for x in S[a_list[0]].mesh]
             N_om = len(Msh)
             DOSproj = [numpy.zeros([ self.shells[ish][3], self.shells[ish][3], N_om ],numpy.float_) for ish in xrange(self.N_shells) ]
             DOS = numpy.zeros([ N_om ],numpy.float_)
             
             ikarray=numpy.array(range(self.Nk))
             stmp = self.add_DC()

             for ik in MPI.slice_array(ikarray):

                 GFsize = [ gf.N1 for sig,gf in S] 
                 unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                       for ib in range(self.NspinblocsGF[self.SO]) ] )

                 if (not unchangedsize):
                     BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                     GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                     a_list = [a for a,al in GFStruct]                                 
                     glist = lambda : [ GFBloc_ReFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]    
                     S = GF(NameList = a_list, BlockList = glist(),Copy=False)
                     S.zero()
                     mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                     for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu
 

                 S <<= GF_Initializers.A_Omega_Plus_B(A=1,B=1j*broadening)
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
                 S *= self.BZ_weights[ik]

                 # non-projected DOS
                 for iom in range(N_om): 
                     for sig,gf in S: DOS[iom] += gf._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                    
                 #projected DOS:
                 if (self.SO==0):
                     for ish in xrange(self.N_shells):
                         tmp = Gproj[ish].copy()
                         for ir in xrange(self.N_parproj[ish]):
                             for sig,gf in tmp: tmp[sig] <<= self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                             Gproj[ish] += tmp
                        
        # collect data from MPI:
        DOS = MPI.all_reduce(MPI.world,DOS,lambda x,y : x+y)
        if (self.SO==0):
            for ish in xrange(self.N_shells):
                Gproj[ish] <<= MPI.all_reduce(MPI.world,Gproj[ish],lambda x,y : x+y)
        MPI.barrier()        

        if (self.SO==0):
            # Symmetrisation:            
            if (self.symm_op!=0): Gproj = self.Symm_par.symmetrise(Gproj)

            # rotation to local coord. system:
            if (self.use_rotations):
                for ish in xrange(self.N_shells):
                    for sig,gf in Gproj[ish]: Gproj[ish][sig] <<= self.rotloc_all(ish,gf,direction='toLocal')
                
            for ish in range(self.N_shells):
                for sig,gf in Gproj[ish]:  DOSproj[ish][:,:,:] += gf._data.array[:,:,:].imag / (-3.1415926535)
	    

        if (MPI.IS_MASTER_NODE()):
            # output to file
            f=open('./DOScorr.dat', 'w')
            for i in range(N_om): f.write("%s    %s\n"%(Msh[i],DOS[i]))
            f.close()    

            if (self.SO==0):
                # partial
                for ish in range(self.N_shells):
                    for i in range(self.shells[ish][3]):
                        for j in range(i,self.shells[ish][3]):
                            Fname = './DOScorr_proj_'+str(ish)+'_'+str(i)+'_'+str(j)+'.dat'
                            f=open(Fname,'w')
                            for iom in range(N_om): f.write("%s    %s\n"%(Msh[iom],DOSproj[ish][i,j,iom]))
                            f.close()




    def calc_FermiSur(self,broadening,Filename, invertAkw=False):
        """ Calculates the FermiSurface with a given Self energy.
            TO BE DONE!!!
            """
        pass        


    def spaghettis(self,broadening,Filename,shift=0.0,plotrange=None, ishell=None, invertAkw=False, Fermisurface=False):
        """ Calculates the correlated band structure with a real-frequency self energy. 
            Filname ... Output file of Wien2K.
            ATTENTION: Many things from the input file are are overwritten!!!"""

        assert type(Filename)==StringType,"Give a valid filename for reading in spaghettis."
        assert hasattr(self,"Sigmaimp"), "Set Sigma First!!"

        assert ((self.SP==0) and (self.SO==0)), "For SP and SO implementation of spaghettis has to be changed!"
        if Fermisurface: ishell=None

        R = Read_Fortran_File2(Filename)
        
        self.Nk = int(R.next())                                 # Number of k-points
        self.N_Orbitals = [ int(R.next()) for i in xrange(self.Nk) ]       # Number of orbitals inside the window 

        
        # Initialise P here for the special k-path:
        self.Proj_Mat = [ [numpy.zeros([self.corr_shells[icrsh][3], self.N_Orbitals[ik]], numpy.complex_) for icrsh in range (self.N_corr_shells)] 
                          for ik in range(self.Nk) ]
 
        for ik in xrange(self.Nk):
            for icrsh in range(self.N_corr_shells):
                                        
                for i in xrange(self.corr_shells[icrsh][3]):    # read real part:
                    for j in xrange(self.N_Orbitals[ik]):
                        self.Proj_Mat[ik][icrsh][i,j] = R.next()
                                
                for i in xrange(self.corr_shells[icrsh][3]):    # read imaginary part:
                    for j in xrange(self.N_Orbitals[ik]):
                        self.Proj_Mat[ik][icrsh][i,j] += 1j * R.next()


        # Hamiltonian for this path:
        self.Hopping = [ numpy.zeros([self.N_Orbitals[ik],self.N_Orbitals[ik]],numpy.complex_) for ik in xrange(self.Nk) ]
        for ik in xrange(self.Nk) :
            for i in xrange(self.N_Orbitals[ik]) :
                for j in xrange(i,self.N_Orbitals[ik]) :
                    self.Hopping[ik][i,j] = R.next() * self.EnergyUnit
		       
            for i in xrange(self.N_Orbitals[ik]) :
                for j in xrange(i,self.N_Orbitals[ik]) :
                    self.Hopping[ik][i,j] += 1j * R.next() * self.EnergyUnit
                    if j !=i : self.Hopping[ik][j,i] = self.Hopping[ik][i,j].conjugate()
        

        # Partial charges projectors:
        self.N_parproj = [int(R.next()) for i in range(self.N_shells)]
        self.Proj_Mat_pc = [ [ [numpy.zeros([self.shells[ish][3], self.N_Orbitals[ik]], numpy.complex_) 
                                for ir in range(self.N_parproj[ish])]
                                for ish in range (self.N_shells)] 
                                for ik in range(self.Nk) ]

        for ish in range(self.N_shells):
            for ik in xrange(self.Nk):
                for ir in range(self.N_parproj[ish]):
                                        
                    for i in xrange(self.shells[ish][3]):    # read real part:
                        for j in xrange(self.N_Orbitals[ik]):
                            self.Proj_Mat_pc[ik][ish][ir][i,j] = R.next()

                            
                    for i in xrange(self.shells[ish][3]):    # read imaginary part:
                        for j in xrange(self.N_Orbitals[ik]):
                            self.Proj_Mat_pc[ik][ish][ir][i,j] += 1j * R.next()


        # Labels:
        self.labels=[]
        # This reads to the end of the file:
        reading = True
        while (reading):
            try:
                R.next() # skip index
                self.labels.append([int(R.next()),R.next()])
            except StopIteration:
                reading=False

        R.close()
        # reading input done!

        if not Fermisurface:
            # print hamiltonian for checks:
            f=open('ham.dat','w')
            for ik in xrange(self.Nk):
                for i in xrange(self.N_Orbitals[ik]):
                    f.write('%s    %s\n'%(ik,self.Hopping[ik][i,i].real))
                f.write('\n')
            f.close()
        
        #=========================================
        # calculate A(k,w):

        if not (ishell is None):
            GFStruct_proj =  [ ('up', range(self.shells[ishell][3])), ('down', range(self.shells[ishell][3])) ] 

        mu = self.Chemical_Potential
        print "Chemical Potential = ",mu

        if self.k_dep_projection == 0:
            assert 0,"Not implemented!!"
        else:
           

            stmp = self.add_DC()
            for ik in xrange(self.Nk):

                if (ik>0):
                    if (self.N_Orbitals[ik]!=self.N_Orbitals[ik-1]):
                        BS = range(self.N_Orbitals[ik])
                        # construct the upfolded Sigma, if #bands changed
                        S = GF(Name_Block_Generator = [ (s,GFBloc_ReFreq(Indices = BS, Mesh = self.Sigmaimp[0].mesh)) for s in ['up','down'] ],Copy = False)
                        mupat = numpy.identity(self.N_Orbitals[ik],numpy.complex_)                                     # change size of mupat
                        mupat *= mu

                else:
                    BS = range(self.N_Orbitals[ik])
                    # construct the upfolded Sigma
                    S = GF(Name_Block_Generator = [ (s,GFBloc_ReFreq(Indices = BS, Mesh = self.Sigmaimp[0].mesh)) for s in ['up','down'] ],Copy = False)
                    mupat = numpy.identity(self.N_Orbitals[ik],numpy.complex_)                                         # construct mupat 
                    mupat *= mu
                    if not (ishell is None):
                        Gproj = GF(Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = S.mesh)) for a,al in GFStruct_proj ], Copy = False)
                        Gproj.zero()
                        
                    # init DOS:
                    M = [x for x in S['up'].mesh]
                    N_om = len(M)
                    if (ishell is None):
                        Akw = numpy.zeros([self.Nk, N_om ],numpy.float_)
                    else:
                        Akw = numpy.zeros([self.shells[ishell][3],self.Nk, N_om ],numpy.float_)
                    if plotrange is None:
                        omminplot = M[0]-0.001
                        ommaxplot = M[N_om-1] + 0.001
                    else:
                        omminplot = plotrange[0]
                        ommaxplot = plotrange[1]
                    if Fermisurface:
                        omminplot = -2.0*broadening
                        ommaxplot =  2.0*broadening
                        Akw = numpy.zeros([self.Nk,1],numpy.float_)
                    

                S <<= GF_Initializers.A_Omega_Plus_B(A=1,B=1j*broadening)
                S -= S.NBlocks * [ self.Hopping[ik] - mupat ]  

                tmp = S.copy()    # init temporary storage
                for icrsh in xrange(self.N_corr_shells):
                    for sig,gf in tmp: 
                        tmp[sig].from_L_G_R(self.Proj_Mat[ik][icrsh].conjugate().transpose(),stmp[icrsh][sig],self.Proj_Mat[ik][icrsh])   
                    S -= tmp      # adding to the upfolded Sigma
                    
                S.invert()
                
                if (ishell is None):
                    # non-projected A(k,w)
                    for iom in range(N_om): 
                        if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                            if Fermisurface:
                                Akw[ik,0] += S['up']._data.array[:,:,iom].imag.trace()/(-3.1415926535) * (M[1]-M[0])
                                Akw[ik,0] += S['down']._data.array[:,:,iom].imag.trace()/(-3.1415926535) * (M[1]-M[0])
                            else:
                                Akw[ik,iom]  = S['up']._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                                Akw[ik,iom] += S['down']._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                                Akw[ik,iom] += ik*shift                       # shift Akw for plotting in xmgrace
                  

                else:
                    # projected A(k,w):
                    Gproj.zero()
                    tmp = Gproj.copy()
                    for ir in xrange(self.N_parproj[ishell]):
                        for sig,gf in tmp: 
                            tmp[sig].from_L_G_R(self.Proj_Mat_pc[ik][ishell][ir],S[sig],self.Proj_Mat_pc[ik][ishell][ir].conjugate().transpose())  # projecting G
                        Gproj += tmp
                    #for sig,gf in Gproj:
                    #    Gproj[sig].from_L_G_R(self.Proj_Mat[ik][0],S[sig],self.Proj_Mat[ik][0].conjugate().transpose())

                    # rotate to local frame
                    if (self.use_rotations):
                        for sig,gf in Gproj: Gproj[sig] <<= self.rotloc(0,gf,direction='toLocal')

                    for iom in range(N_om): 
                        if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                            for ish in range(self.shells[ishell][3]):
                                Akw[ish,ik,iom]  = Gproj['up']._data.array[ish,ish,iom].imag/(-3.1415926535)
                                Akw[ish,ik,iom] += Gproj['down']._data.array[ish,ish,iom].imag/(-3.1415926535)
                   
            
            # END k-LOOP

        if (ishell is None):
                            
            if (invertAkw):
                maxAkw=Akw.max()
                minAkw=Akw.min()
                
            # open file for storage:
            if Fermisurface:
                f=open('FS.dat','w')
            else:
                f=open('Akw.dat','w')

            for ik in range(self.Nk):
                if Fermisurface:
                    if (invertAkw):
                        Akw[ik,0] = 1.0/(minAkw-maxAkw)*(Akw[ik,0] - maxAkw)
                    f.write('%s    %s\n'%(ik,Akw[ik,0]))
                else:
                    for iom in range(N_om): 
                        if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                            if (invertAkw):
                                Akw[ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ik,iom] - maxAkw)
                            if (shift>0.0001):
                                f.write('%s      %s\n'%(M[iom],Akw[ik,iom]))
                            else:
                                f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ik,iom]))

                    f.write('\n')
 
            f.close()

        else:

            for ish in range(self.shells[ishell][3]):
                    
                if (invertAkw):
                    maxAkw=Akw[ish,:,:].max()
                    minAkw=Akw[ish,:,:].min()

                f=open('Akw_proj'+str(ish)+'.dat','w') 

                for ik in range(self.Nk):
                    for iom in range(N_om): 
                        if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                            if (invertAkw):
                                Akw[ish,ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ish,ik,iom] - maxAkw)
                            if (shift>0.0001):
                                f.write('%s      %s\n'%(M[iom],Akw[ish,ik,iom]))
                            else:
                                f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ish,ik,iom]))

                    f.write('\n')
 
                f.close()



    def corr_bands(self,broadening,Filename,Z=None,Deltamu=None,shift=0.0,plotrange=None,invertAkw=True, orbtoproj = False):

        assert type(Filename)==StringType,"Give a valid filename for reading in spaghettis."
        

        R = Read_Fortran_File2(Filename)
        
        self.Nk = int(R.next())                                 # Number of k-points
        self.N_Orbitals = [ int(R.next()) for i in xrange(self.Nk) ]       # Number of orbitals inside the window 

        
        # Initialise P here for the special k-path:
        self.Proj_Mat = [ [numpy.zeros([self.corr_shells[icrsh][3], self.N_Orbitals[ik]], numpy.complex_) for icrsh in range (self.N_corr_shells)] 
                          for ik in range(self.Nk) ]
 
        for ik in xrange(self.Nk):
            for icrsh in range(self.N_corr_shells):
                                        
                for i in xrange(self.corr_shells[icrsh][3]):    # read real part:
                    for j in xrange(self.N_Orbitals[ik]):
                        self.Proj_Mat[ik][icrsh][i,j] = R.next()
                                
                for i in xrange(self.corr_shells[icrsh][3]):    # read imaginary part:
                    for j in xrange(self.N_Orbitals[ik]):
                        self.Proj_Mat[ik][icrsh][i,j] += 1j * R.next()


        # Hamiltonian for this path:
        self.Hopping = [ numpy.zeros([self.N_Orbitals[ik],self.N_Orbitals[ik]],numpy.complex_) for ik in xrange(self.Nk) ]
        for ik in xrange(self.Nk) :
            for i in xrange(self.N_Orbitals[ik]) :
                for j in xrange(i,self.N_Orbitals[ik]) :
                    self.Hopping[ik][i,j] = R.next() * self.EnergyUnit
		       
            for i in xrange(self.N_Orbitals[ik]) :
                for j in xrange(i,self.N_Orbitals[ik]) :
                    self.Hopping[ik][i,j] += 1j * R.next() * self.EnergyUnit
                    if j !=i : self.Hopping[ik][j,i] = self.Hopping[ik][i,j].conjugate()
        

        # Partial charges projectors:
        self.N_parproj = [int(R.next()) for i in range(self.N_shells)]
        self.Proj_Mat_pc = [ [ [numpy.zeros([self.shells[ish][3], self.N_Orbitals[ik]], numpy.complex_) 
                                for ir in range(self.N_parproj[ish])]
                                for ish in range (self.N_shells)] 
                                for ik in range(self.Nk) ]

        for ish in range(self.N_shells):
            for ik in xrange(self.Nk):
                for ir in range(self.N_parproj[ish]):
                                        
                    for i in xrange(self.shells[ish][3]):    # read real part:
                        for j in xrange(self.N_Orbitals[ik]):
                            self.Proj_Mat_pc[ik][ish][ir][i,j] = R.next()

                            
                    for i in xrange(self.shells[ish][3]):    # read imaginary part:
                        for j in xrange(self.N_Orbitals[ik]):
                            self.Proj_Mat_pc[ik][ish][ir][i,j] += 1j * R.next()


        # Labels:
        self.labels=[]
        # This reads to the end of the file:
        reading = True
        while (reading):
            try:
                R.next() # skip index
                self.labels.append([int(R.next()),R.next()])
            except StopIteration:
                reading=False

        R.close()
        # reading input done!

        # print hamiltonian for checks:
        f=open('ham.dat','w')
        for ik in xrange(self.Nk):
            for i in xrange(self.N_Orbitals[ik]):
                f.write('%s    %s\n'%(ik,self.Hopping[ik][i,i].real))
            f.write('\n')
        f.close()

        

        mu = self.Chemical_Potential
            
           
        for ik in xrange(self.Nk):

            if (ik>0):
                if (self.N_Orbitals[ik]!=self.N_Orbitals[ik-1]):
                    BS = range(self.N_Orbitals[ik])
                    # construct the upfolded Sigma, if #bands changed
                    S = GF(Name_Block_Generator = [ (s,GFBloc_ReFreq(Indices = BS, Mesh = self.Sigmaimp[0].mesh)) for s in ['up','down'] ],Copy = False)
                    mupat = numpy.identity(self.N_Orbitals[ik],numpy.complex_)                                     # change size of mupat
                    mupat *= mu
                    
            else:
                BS = range(self.N_Orbitals[ik])
                # construct the upfolded Sigma
                S = GF(Name_Block_Generator = [ (s,GFBloc_ReFreq(Indices = BS, Mesh = self.Sigmaimp[0].mesh)) for s in ['up','down'] ],Copy = False)
                mupat = numpy.identity(self.N_Orbitals[ik],numpy.complex_)                                         # construct mupat 
                mupat *= mu

                Gproj = [ self.Sigmaimp[icrsh].copy() for icrsh in xrange(self.N_corr_shells) ]   # this list will be returned  
                for icrsh in range(self.N_corr_shells): Gproj[icrsh].zero()
                # init DOS:
                M = [x for x in S['up'].mesh]
                N_om = len(M)
                if (orbtoproj == False):
                    Akw = numpy.zeros([self.Nk, N_om ],numpy.float_)
                else:
                    Akw = numpy.zeros([self.corr_shells[0][3],self.Nk,N_om], numpy.float_)

                if plotrange is None:
                    omminplot = M[0]-0.001
                    ommaxplot = M[N_om-1] + 0.001
                else:
                    omminplot = plotrange[0]
                    ommaxplot = plotrange[1]
                        
                       

            if (Z is None):
                A = 1.0
            else:
                # set the quasiparticle renormalization:
                A = numpy.zeros([self.N_Orbitals[ik],self.N_Orbitals[ik]], numpy.float_)
                for icrsh in range(self.N_corr_shells):
                    A += numpy.dot(self.Proj_Mat[ik][icrsh].conjugate().transpose(),numpy.dot(Z,self.Proj_Mat[ik][icrsh]))
                for nu in range(self.N_Orbitals[ik]): A[nu,nu] += 1.0

            Dmu = numpy.zeros([self.N_Orbitals[ik],self.N_Orbitals[ik]], numpy.float_)
            if not (Deltamu is None):
                for icrsh in range(self.N_corr_shells):
                    Dmu += numpy.dot(self.Proj_Mat[ik][icrsh].conjugate().transpose(),numpy.dot(Deltamu,self.Proj_Mat[ik][icrsh]))
            
                
            S <<= GF_Initializers.A_Omega_Plus_B(A=A,B=1j*broadening)
            S -= S.NBlocks * [ self.Hopping[ik] - mupat + Dmu ]  


                    
            S.invert()
                
            # projected A(k,w):
            for icrsh in range(self.N_corr_shells):
                for sig,gf in Gproj[icrsh]:
                    Gproj[icrsh][sig].from_L_G_R(self.Proj_Mat[ik][icrsh],S[sig],self.Proj_Mat[ik][icrsh].conjugate().transpose())
              
            
            # rotate to local frame:
            if (self.use_rotations):
                for icrsh in xrange(self.N_corr_shells):
                    for sig,gf in Gproj[icrsh]: Gproj[icrsh][sig] <<= self.rotloc(icrsh,gf,direction='toLocal')
            
            
            for iom in range(N_om): 
                if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                    if (orbtoproj==False):
                        Akw[ik,iom]  = Gproj[0]['up']._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                        Akw[ik,iom] += Gproj[0]['down']._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                    else:
                        for ish in range(self.corr_shells[0][3]):
                            Akw[ish,ik,iom]  = Gproj[0]['up']._data.array[ish,ish,iom].imag/(-3.1415926535)
                            Akw[ish,ik,iom] += Gproj[0]['down']._data.array[ish,ish,iom].imag/(-3.1415926535)
 
        # END k-LOOP

        if (orbtoproj==False):
            if (invertAkw):
                maxAkw=Akw.max()
                minAkw=Akw.min()
                
            # open file for storage:
            f=open('Akw_corrbands.dat','w')

            for ik in range(self.Nk):
                for iom in range(N_om): 
                    if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                        if (invertAkw):
                            Akw[ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ik,iom] - maxAkw)
                        if (shift>0.0001):
                            f.write('%s      %s\n'%(M[iom],Akw[ik,iom]))
                        else:
                            f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ik,iom]))

                f.write('\n')
 
            f.close()

        else:

            for ish in range(self.corr_shells[0][3]):

                if (invertAkw):
                    maxAkw=Akw[ish,:,:].max()
                    minAkw=Akw[ish,:,:].min()

                f=open('Akw_corrbands_proj'+str(ish)+'.dat','w')

                for ik in range(self.Nk):
                    for iom in range(N_om): 
                        if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                            if (invertAkw):
                                Akw[ish,ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ish,ik,iom] - maxAkw)
                            if (shift>0.0001):
                                f.write('%s      %s\n'%(M[iom],Akw[ish,ik,iom]))
                            else:
                                f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ish,ik,iom]))

                    f.write('\n')
 
                f.close()
                


    
    def constr_Sigma_ME(self,Filename, Beta, N_om, orb = 0):
        """Uses Data from files to construct a GF object on the real axis."""

        
        #first get the mesh out of one of the files:
        if (len(self.GFStruct_Solver[orb][0][1])==1):
            Fname = Filename+'_'+self.GFStruct_Solver[orb][0][0]+'.dat'
        else:
            Fname = Filename+'_'+self.GFStruct_Solver[orb][0][0]+'/'+str(self.GFStruct_Solver[orb][0][1][0])+'_'+str(self.GFStruct_Solver[orb][0][1][0])+'.dat'

        R = Read_Fortran_File(Fname)
        mesh = numpy.zeros([N_om],numpy.float_)
        try:
            for i in xrange(N_om): 
                mesh[i] = R.next()
                sk = R.next()
                sk = R.next()
                
        except StopIteration : # a more explicit error if the file is corrupted.
            raise "SumK_LDA.read_Sigma_ME : reading file failed!"
        R.close()

        # now initialize the GF with the mesh
        a_list = [a for a,al in self.GFStruct_Solver[orb]]
        glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = mesh) for a,al in self.GFStruct_Solver[orb] ] 
        SigmaME = GF(NameList = a_list, BlockList = glist(),Copy=False)
        SigmaME.load(Filename)
        SigmaME.Note='ReFreq'          # This is important for the put_Sigma routine!!!

        return SigmaME


        

    def partial_charges(self,Proj_filename=False,Symm_file='Symmetry.dat'):
        """calculates the density matrix in the desired shells. The projectors to these shells are given in the file.
           It is assumed that the N_Orbitals is the same as read in __init__ """

        
        if (Proj_filename):
            # if a filename is given, read the Projectors and the Symmetry operations
            self.__readprojfile(Proj_filename = Proj_filename, Symm_file = Symm_file)
            
        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]
        # Desnity matrix in the window
        self.Dens_Mat_window = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.N_shells)]
                                 for isp in range(len(bln)) ]  # arrays for all blocs

        if hasattr(self,"Sigmaimp"):
            self.Dens_Mat_window = self.__pcwithSigma()
        else:
            self.Dens_Mat_window = self.__pcwoSigma()

            
        
        # add Density matrices to get the total:
        Dens_Mat = [ [ self.Dens_Mat_below[ntoi[bln[isp]]][ish]+self.Dens_Mat_window[isp][ish] for ish in range(self.N_shells)]
                     for isp in range(len(bln)) ]
        

        return Dens_Mat


    def __readprojfile(self,Proj_filename,Symm_file):
        
        # Density matrix below the energy window, for information:
        self.Dens_Mat_below = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.N_shells)] 
                                for isp in range(self.Nspinblocs) ]

        R = Read_Fortran_File(Proj_filename)
        try:

            if self.k_dep_projection == 0:
                # Projection Matrix does not depend on k:
                # i.e., number of bands in the energy window is the same for all k!!!
                # read it once and for all

                assert 0, "Partial charges only implemented for k-dependent projectors!!!"
                   
            else:
                
                # read the number of projectors for each shell:
                self.N_parproj = [int(R.next()) for i in range(self.N_shells)]
                
                    # Initialise P, here a double list of matrices:
                self.Proj_Mat_pc = [ [ [ [numpy.zeros([self.shells[ish][3], self.N_Orbitals[ik][isp]], numpy.complex_) 
                                          for ir in range(self.N_parproj[ish])]
                                         for ish in range (self.N_shells) ]
                                       for isp in range(self.Nspinblocs) ]
                                     for ik in range(self.Nk) ]
                                     
 
                # Global -> local rotation matrices:
                if (self.use_rotations):
                    self.rotmat_all = [numpy.identity(self.shells[ish][3],numpy.complex_) for ish in xrange(self.N_shells)]

                for isp in range(self.Nspinblocs):
                    for ish in range(self.N_shells):

                        # read first the projectors for this orbital:
                        for ik in xrange(self.Nk):
                            for ir in range(self.N_parproj[ish]):
                                        
                                for i in xrange(self.shells[ish][3]):    # read real part:
                                    for j in xrange(self.N_Orbitals[ik][isp]):
                                        self.Proj_Mat_pc[ik][isp][ish][ir][i,j] = R.next()

                            
                                for i in xrange(self.shells[ish][3]):    # read imaginary part:
                                    for j in xrange(self.N_Orbitals[ik][isp]):
                                        self.Proj_Mat_pc[ik][isp][ish][ir][i,j] += 1j * R.next()
                                        
                        # now read the Density Matrix for this orbital below the energy window:
                        for i in xrange(self.shells[ish][3]):    # read real part:
                            for j in xrange(self.shells[ish][3]):
                                self.Dens_Mat_below[isp][ish][j,i] = R.next()
                        for i in xrange(self.shells[ish][3]):    # read imaginary part:
                            for j in xrange(self.shells[ish][3]):
                                self.Dens_Mat_below[isp][ish][j,i] += 1j * R.next()
                        if (self.Nspinblocs==1): self.Dens_Mat_below[isp][ish] /= 2.0

                        # Global -> local rotation matrix for this shell:
                        if (self.use_rotations):
                            if (self.shells[ish][2]==0): # read p, d, and f matrices, s is unity
                                self.rotmat_all[ish][0,0] = 1.0
                            else:
                                for i in xrange(self.shells[ish][3]):    # read real part:
                                    for j in xrange(self.shells[ish][3]):
                                        self.rotmat_all[ish][j,i] = R.next()
                                for i in xrange(self.shells[ish][3]):    # read imaginary part:
                                    for j in xrange(self.shells[ish][3]):
                                        self.rotmat_all[ish][j,i] += 1j * R.next()


        except StopIteration : # a more explicit error if the file is corrupted.
            raise "SumK_LDA.partial_charges : reading file for Projectors failed!"

        R.close()
        self.__Proj_Mat_pc_read = True

        # initialize the symmetrisation routine for partial charges:
        if (self.symm_op!=0): self.Symm_par = Symmetry(Symm_file,self.shells)



    def __pcwithSigma(self):
        """calculates the density matrix, if Sigmaimp has been set"""

        #print "use Sigma\n"

        mu = self.Chemical_Potential
        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]
        dm = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.N_shells)]   
               for isp in range(len(bln)) ]    # init the density matrix

        GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in bln ]  for i in xrange(self.N_shells) ]
        Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct_proj[ish] ], Copy = False)
                 for ish in xrange(self.N_shells)]
        for ish in xrange(self.N_shells): Gproj[ish].zero()

        BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
        GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
        a_list = [a for a,al in GFStruct]  
        glist = lambda : [ GFBloc_ImFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]  
        S = GF(NameList = a_list, BlockList = glist(),Copy=False)
        S.zero()

        mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln]   # construct mupat
        for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu

        ikarray=numpy.array(range(self.Nk))
        stmp = self.add_DC()
        for ik in MPI.slice_array(ikarray):
            print ik

            GFsize = [ gf.N1 for sig,gf in S] 
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                       for ib in range(self.NspinblocsGF[self.SO]) ] )

            if (not unchangedsize):
                BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                a_list = [a for a,al in GFStruct]                                 
                glist = lambda : [ GFBloc_ImFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]    
                S = GF(NameList = a_list, BlockList = glist(),Copy=False)
                S.zero()
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
            S *= self.BZ_weights[ik]

            for ish in xrange(self.N_shells):
                tmp = Gproj[ish].copy()
                for ir in xrange(self.N_parproj[ish]):
                    for sig,gf in tmp: tmp[sig] <<= self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                    Gproj[ish] += tmp

        #collect data from MPI:
        for ish in xrange(self.N_shells):
            Gproj[ish] <<= MPI.all_reduce(MPI.world,Gproj[ish],lambda x,y : x+y)
        MPI.barrier()

        # Symmetrisation:
        if (self.symm_op!=0): Gproj = self.Symm_par.symmetrise(Gproj)
          
        
        for ish in xrange(self.N_shells):
                
            # Rotation to local:
            if (self.use_rotations):
                for sig,gf in Gproj[ish]: Gproj[ish][sig] <<= self.rotloc_all(ish,gf,direction='toLocal')
                
            isp = 0
            for sig,gf in Gproj[ish]: #dmg.append(Gproj[ish].density()[sig])
                dm[isp][ish] = Gproj[ish].density()[sig]
                isp+=1
            


        return dm


    def __pcwoSigma(self):
        """calculates the density matrix, if Sigmaimp has been set"""

        #print "use Sigma\n"
        
        Beta = 40.0

        mu = self.Chemical_Potential
        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]
        dm = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.N_shells)]   
               for isp in range(len(bln)) ]    # init the density matrix

        BS = [ range(self.N_Orbitals[0][ntoi[ib]]) for ib in bln ]
        GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
        a_list = [a for a,al in GFStruct]  
        glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta=Beta) for a,al in GFStruct]  
        S = GF(NameList = a_list, BlockList = glist(),Copy=False)
        S.zero()
        mupat = [numpy.identity(self.N_Orbitals[0][ntoi[bl]],numpy.complex_) for bl in bln]   # construct mupat
        for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu

        GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in bln ]  for i in xrange(self.N_shells) ]
        Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = S.mesh)) for a,al in GFStruct_proj[ish] ], Copy = False)
                 for ish in xrange(self.N_shells)]
        for ish in xrange(self.N_shells): Gproj[ish].zero()


        ikarray=numpy.array(range(self.Nk))
        
        for ik in MPI.slice_array(ikarray):
       
            GFsize = [ gf.N1 for sig,gf in S] 
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                                       for ib in range(self.NspinblocsGF[self.SO]) ] )

            if (not unchangedsize):
                BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
                GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
                a_list = [a for a,al in GFStruct]                                 
                glist = lambda : [ GFBloc_ImFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]    
                S = GF(NameList = a_list, BlockList = glist(),Copy=False)
                S.zero()
                mupat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]   # change size of mupat
                for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu

            S <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
            M = copy.deepcopy(mupat)
            for ibl in range(self.NspinblocsGF[self.SO]): 
                ind = ntoi[bln[ibl]]
                M[ibl] = self.Hopping[ik][ind] - mupat[ibl]
            S -= M  

            S.invert()
            S *= self.BZ_weights[ik]

            for ish in xrange(self.N_shells):
                tmp = Gproj[ish].copy()
                for ir in xrange(self.N_parproj[ish]):
                    for sig,gf in tmp: tmp[sig] <<= self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                    Gproj[ish] += tmp

        #collect data from MPI:
        for ish in xrange(self.N_shells):
            Gproj[ish] <<= MPI.all_reduce(MPI.world,Gproj[ish],lambda x,y : x+y)
        MPI.barrier()

        # Symmetrisation:
        if (self.symm_op!=0): Gproj = self.Symm_par.symmetrise(Gproj)
          
        
        for ish in xrange(self.N_shells):
                
            # Rotation to local:
            if (self.use_rotations):
                for sig,gf in Gproj[ish]: Gproj[ish][sig] <<= self.rotloc_all(ish,gf,direction='toLocal')
                
            isp = 0
            for sig,gf in Gproj[ish]: #dmg.append(Gproj[ish].density()[sig])
                dm[isp][ish] = Gproj[ish].density()[sig]
                isp+=1
            


        return dm



