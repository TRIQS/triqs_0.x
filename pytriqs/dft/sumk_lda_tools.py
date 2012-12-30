
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
import pytriqs.base.utility.Dichotomy as Dichotomy
from pytriqs.base.gf_local.block_gf import GF
from pytriqs.base.gf_local.gf_imfreq import GFBloc_ImFreq
from pytriqs.base.gf_local.gf_refreq import GFBloc_ReFreq
from pytriqs.base.gf_local.gf_imtime import GFBloc_ImTime
from pytriqs.base.gf_local import gf_init
from pytriqs.solvers.operators import *
from pytriqs.base.utility.myUtils import Sum
import pytriqs.base.utility.MPI as MPI
from datetime import datetime

from pytriqs.dft.symmetry import *
from pytriqs.dft.sumk_lda import SumK_LDA

import string


def Read_Fortran_File (filename):
    """ Returns a generator that yields all numbers in the Fortran file as float, one by one"""
    import os.path
    if not(os.path.exists(filename)) : raise IOError, "File %s does not exists"%filename
    for line in open(filename,'r') :
	for x in line.replace('D','E').split() : 
	    yield string.atof(x)


class SumK_LDA_tools(SumK_LDA):
    """Extends the SumK_LDA class with some tools for analysing the data."""


    def __init__(self, HDFfile, mu = 0.0, hfield = 0.0, UseLDABlocs = False, LDAdata = 'SumK_LDA', Symmcorrdata = 'SymmCorr',
                 ParProjdata = 'SumK_LDA_ParProj', Symmpardata = 'SymmPar', Bandsdata = 'SumK_LDA_Bands'):

        self.Gupf_refreq = None
        SumK_LDA.__init__(self,HDFfile=HDFfile,mu=mu,hfield=hfield,UseLDABlocs=UseLDABlocs,LDAdata=LDAdata,
                          Symmcorrdata=Symmcorrdata,ParProjdata=ParProjdata,Symmpardata=Symmpardata,
                          Bandsdata=Bandsdata)   
    

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

        #gf_rotated = gf_to_rotate.copy()
        #if (direction=='toGlobal'):
        #    gf_rotated.from_L_G_R(self.rotmat_all[ish],gf_to_rotate,self.rotmat_all[ish].conjugate().transpose())
        #elif (direction=='toLocal'):
        #    gf_rotated.from_L_G_R(self.rotmat_all[ish].conjugate().transpose(),gf_to_rotate,self.rotmat_all[ish])

        gf_rotated = gf_to_rotate.copy()
        if (direction=='toGlobal'):
            #if (self.rotmat_timeinv[ish]==1): gf_rotated <<= gf_rotated.transpose()
            #gf_rotated.from_L_G_R(self.rotmat[ish].transpose(),gf_rotated,self.rotmat[ish].conjugate())
            if ((self.rotmat_all_timeinv[ish]==1) and (self.SO)):
                gf_rotated <<= gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rotmat_all[ish].conjugate(),gf_rotated,self.rotmat_all[ish].transpose())
            else:
                gf_rotated.from_L_G_R(self.rotmat_all[ish],gf_rotated,self.rotmat_all[ish].conjugate().transpose())
            
        elif (direction=='toLocal'):
            if ((self.rotmat_all_timeinv[ish]==1)and(self.SO)):
                gf_rotated <<= gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rotmat_all[ish].transpose(),gf_rotated,self.rotmat_all[ish].conjugate())
            else:
                gf_rotated.from_L_G_R(self.rotmat_all[ish].conjugate().transpose(),gf_rotated,self.rotmat_all[ish])
                
    
        return gf_rotated


    def latticeGF_realfreq(self, ik, mu, broadening, Mesh=None, Beta=40, withSigma=True):
        """Calculates the lattice Green function on the real frequency axis. If self energy is
           present and withSigma=True, the mesh is taken from Sigma. Otherwise, the mesh has to be given."""

        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        if (not hasattr(self,"Sigmaimp")): withSigma=False
        if (withSigma): 
            assert self.Sigmaimp[0].Note=='ReFreq',"Real frequency Sigma needed for latticeGF_realfreq!"
            Beta = self.Sigmaimp[0].Beta
            stmp = self.add_DC()
        else:
            assert (not (Mesh is None)),"Without Sigma, give the mesh for latticeGF_realfreq!"

	if (self.Gupf_refreq is None):
            # first setting up of Gupf_refreq
            BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
            GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
            a_list = [a for a,al in GFStruct]
            if (withSigma):
                glist = lambda : [ GFBloc_ReFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]
            else:
                glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in GFStruct] 
            self.Gupf_refreq = GF(NameList = a_list, BlockList = glist(),Copy=False)
            self.Gupf_refreq.zero()

        GFsize = [ gf.N1 for sig,gf in self.Gupf_refreq]
        unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib]
                               for ib in range(self.NspinblocsGF[self.SO]) ] )

        if (not unchangedsize):
            BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
            GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
            a_list = [a for a,al in GFStruct]
            if (withSigma):
                glist = lambda : [ GFBloc_ReFreq(Indices = al, Mesh=self.Sigmaimp[0].mesh) for a,al in GFStruct]
            else:
                glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in GFStruct]
            self.Gupf_refreq = GF(NameList = a_list, BlockList = glist(),Copy=False)
            self.Gupf_refreq.zero()
        
        idmat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]

        self.Gupf_refreq <<= gf_init.A_Omega_Plus_B(A=1,B=1j*broadening)
        M = copy.deepcopy(idmat)
        for ibl in range(self.NspinblocsGF[self.SO]):
            ind = ntoi[bln[ibl]]
            M[ibl] = self.Hopping[ik][ind] - (idmat[ibl]*mu) - (idmat[ibl] * self.hfield * (1-2*ibl))
        self.Gupf_refreq -= M

        if (withSigma):
            tmp = self.Gupf_refreq.copy()    # init temporary storage
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in tmp: tmp[sig] <<= self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)
                self.Gupf_refreq -= tmp      # adding to the upfolded GF

        self.Gupf_refreq.invert()

        return self.Gupf_refreq



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
                dl = self.corr_shells[self.invshellmap[icrsh]][3]
                DOSproj[icrsh][bn] = numpy.zeros([N_om],numpy.float_)
                DOSproj_orb[icrsh][bn] = numpy.zeros([dl,dl,N_om],numpy.float_)

      
        for i in range(N_om): Mesh[i] = ommin + delta_om * i

        # init:
        Gloc = []
        for icrsh in range(self.N_corr_shells):
            b_list = [a for a,al in self.GFStruct_corr[icrsh]]
            glist = lambda : [ GFBloc_ReFreq(Indices = al, Beta = Beta, MeshArray = Mesh) for a,al in self.GFStruct_corr[icrsh]]   
            Gloc.append(GF(NameList = b_list, BlockList = glist(),Copy=False))
        for icrsh in xrange(self.N_corr_shells): Gloc[icrsh].zero()                        # initialize to zero
            
        for ik in xrange(self.Nk):

            Gupf=self.latticeGF_realfreq(ik=ik,mu=self.Chemical_Potential,broadening=broadening,Beta=Beta,Mesh=Mesh,withSigma=False)
            Gupf *= self.BZ_weights[ik]

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
     
        # output:
        if (MPI.IS_MASTER_NODE()):
            for bn in self.blocnames[self.SO]:
                f=open('DOS%s.dat'%bn, 'w')
                for i in range(N_om): f.write("%s    %s\n"%(Mesh[i],DOS[bn][i]))
                f.close()  

                for ish in range(self.N_inequiv_corr_shells):
                    f=open('DOS%s_proj%s.dat'%(bn,ish),'w')
                    for i in range(N_om): f.write("%s    %s\n"%(Mesh[i],DOSproj[ish][bn][i]))
                    f.close()  
            
                    for i in range(self.corr_shells[self.invshellmap[ish]][3]):
                        for j in range(i,self.corr_shells[self.invshellmap[ish]][3]):
                            Fname = 'DOS'+bn+'_proj'+str(ish)+'_'+str(i)+'_'+str(j)+'.dat'
                            f=open(Fname,'w')
                            for iom in range(N_om): f.write("%s    %s\n"%(Mesh[iom],DOSproj_orb[ish][bn][i,j,iom]))
                            f.close()
 

       

    def read_ParProj_input_from_HDF(self):
        """
        Reads the data for the partial projectors from the HDF file
        """

        thingstoread = ['Dens_Mat_below','N_parproj','Proj_Mat_pc','rotmat_all','rotmat_all_timeinv']
        retval = self.read_input_from_HDF(SubGrp=self.ParProjdata,thingstoread = thingstoread)
        return retval



    def DOSpartial(self,broadening=0.01):
        """calculates the orbitally-resolved DOS"""

        assert hasattr(self,"Sigmaimp"), "Set Sigma First!!"

        #thingstoread = ['Dens_Mat_below','N_parproj','Proj_Mat_pc','rotmat_all']
        #retval = self.read_input_from_HDF(SubGrp=self.ParProjdata, thingstoread=thingstoread)
        retval = self.read_ParProj_input_from_HDF()
        if not retval: return retval
        if self.symm_op: self.Symm_par = Symmetry(self.HDFfile,subgroup=self.Symmpardata)

        mu = self.Chemical_Potential

        GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in self.blocnames[self.SO] ]  for i in xrange(self.N_shells) ]
        Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct_proj[ish] ], Copy = False ) 
                 for ish in xrange(self.N_shells)]
        for ish in range(self.N_shells): Gproj[ish].zero()

        Msh = [x for x in self.Sigmaimp[0].mesh]
        N_om = len(Msh)

        DOS = {}
        for bn in self.blocnames[self.SO]:
            DOS[bn] = numpy.zeros([N_om],numpy.float_)

        DOSproj     = [ {} for ish in range(self.N_shells) ]
        DOSproj_orb = [ {} for ish in range(self.N_shells) ]
        for ish in range(self.N_shells):
            for bn in self.blocnames[self.SO]:
                dl = self.shells[ish][3]
                DOSproj[ish][bn] = numpy.zeros([N_om],numpy.float_)
                DOSproj_orb[ish][bn] = numpy.zeros([dl,dl,N_om],numpy.float_)

        ikarray=numpy.array(range(self.Nk))

        for ik in MPI.slice_array(ikarray):

            S = self.latticeGF_realfreq(ik=ik,mu=mu,broadening=broadening)
            S *= self.BZ_weights[ik]

            # non-projected DOS
            for iom in range(N_om): 
                for sig,gf in S: DOS[sig][iom] += gf._data.array[:,:,iom].imag.trace()/(-3.1415926535)
               
            #projected DOS:
            for ish in xrange(self.N_shells):
                tmp = Gproj[ish].copy()
                for ir in xrange(self.N_parproj[ish]):
                    for sig,gf in tmp: tmp[sig] <<= self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                    Gproj[ish] += tmp
                   
        # collect data from MPI:
        for sig in DOS:
            DOS[sig] = MPI.all_reduce(MPI.world,DOS[sig],lambda x,y : x+y)
        for ish in xrange(self.N_shells):
            Gproj[ish] <<= MPI.all_reduce(MPI.world,Gproj[ish],lambda x,y : x+y)
        MPI.barrier()        
                  
        if (self.symm_op!=0): Gproj = self.Symm_par.symmetrise(Gproj)

        # rotation to local coord. system:
        if (self.use_rotations):
            for ish in xrange(self.N_shells):
                for sig,gf in Gproj[ish]: Gproj[ish][sig] <<= self.rotloc_all(ish,gf,direction='toLocal')
                
        for ish in range(self.N_shells):
            for sig,gf in Gproj[ish]:  
                for iom in range(N_om): DOSproj[ish][sig][iom] += gf._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                DOSproj_orb[ish][sig][:,:,:] += gf._data.array[:,:,:].imag / (-3.1415926535)
	    

        if (MPI.IS_MASTER_NODE()):
            # output to files
            for bn in self.blocnames[self.SO]:
                f=open('./DOScorr%s.dat'%bn, 'w')
                for i in range(N_om): f.write("%s    %s\n"%(Msh[i],DOS[bn][i]))
                f.close()    

                # partial
                for ish in range(self.N_shells):
                    f=open('DOScorr%s_proj%s.dat'%(bn,ish),'w')
                    for i in range(N_om): f.write("%s    %s\n"%(Msh[i],DOSproj[ish][bn][i]))
                    f.close()
 
                    for i in range(self.shells[ish][3]):
                        for j in range(i,self.shells[ish][3]):
                            Fname = './DOScorr'+bn+'_proj'+str(ish)+'_'+str(i)+'_'+str(j)+'.dat'
                            f=open(Fname,'w')
                            for iom in range(N_om): f.write("%s    %s\n"%(Msh[iom],DOSproj_orb[ish][bn][i,j,iom]))
                            f.close()




    def spaghettis(self,broadening,shift=0.0,plotrange=None, ishell=None, invertAkw=False, Fermisurface=False):
        """ Calculates the correlated band structure with a real-frequency self energy. 
            ATTENTION: Many things from the original input file are are overwritten!!!"""

        assert hasattr(self,"Sigmaimp"), "Set Sigma First!!"
        thingstoread = ['Nk','N_Orbitals','Proj_Mat','Hopping','N_parproj','Proj_Mat_pc']
        retval = self.read_input_from_HDF(SubGrp=self.Bandsdata,thingstoread=thingstoread)
        if not retval: return retval

        if Fermisurface: ishell=None
       
        # print hamiltonian for checks:
        if ((self.SP==1)and(self.SO==0)):
            f1=open('hamup.dat','w')
            f2=open('hamdn.dat','w')
           
            for ik in xrange(self.Nk): 
                for i in xrange(self.N_Orbitals[ik][0]):
                    f1.write('%s    %s\n'%(ik,self.Hopping[ik][0][i,i].real))
                for i in xrange(self.N_Orbitals[ik][1]):
                    f2.write('%s    %s\n'%(ik,self.Hopping[ik][1][i,i].real))
                f1.write('\n')
                f2.write('\n')
            f1.close()
            f2.close()
        else:
            f=open('ham.dat','w')
            for ik in xrange(self.Nk):
                for i in xrange(self.N_Orbitals[ik][0]):
                    f.write('%s    %s\n'%(ik,self.Hopping[ik][0][i,i].real))
                f.write('\n')
            f.close()

        
        #=========================================
        # calculate A(k,w):

        mu = self.Chemical_Potential
        bln = self.blocnames[self.SO]

        # init DOS:
        M = [x for x in self.Sigmaimp[0].mesh]
        N_om = len(M)

        if plotrange is None:
            omminplot = M[0]-0.001
            ommaxplot = M[N_om-1] + 0.001
        else:
            omminplot = plotrange[0]
            ommaxplot = plotrange[1]

        if (ishell is None):
            Akw = {}
            for ibn in bln: Akw[ibn] = numpy.zeros([self.Nk, N_om ],numpy.float_)
        else:
            Akw = {}
            for ibn in bln: Akw[ibn] = numpy.zeros([self.shells[ishell][3],self.Nk, N_om ],numpy.float_)

        if Fermisurface:
            omminplot = -2.0*broadening
            ommaxplot =  2.0*broadening
            Akw = {}
            for ibn in bln: Akw[ibn] = numpy.zeros([self.Nk,1],numpy.float_)

        if not (ishell is None):
            GFStruct_proj =  [ (al, range(self.shells[ishell][3])) for al in bln ]
            Gproj = GF(Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct_proj ], Copy = False)
            Gproj.zero()

        for ik in xrange(self.Nk):

            S = self.latticeGF_realfreq(ik=ik,mu=mu,broadening=broadening)               
            if (ishell is None):
                # non-projected A(k,w)
                for iom in range(N_om): 
                    if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                        if Fermisurface:
                            for sig,gf in S: Akw[sig][ik,0] += gf._data.array[:,:,iom].imag.trace()/(-3.1415926535) * (M[1]-M[0])
                        else:
                            for sig,gf in S: Akw[sig][ik,iom] += gf._data.array[:,:,iom].imag.trace()/(-3.1415926535)
                            Akw[sig][ik,iom] += ik*shift                       # shift Akw for plotting in xmgrace
                  

            else:
                # projected A(k,w):
                Gproj.zero()
                tmp = Gproj.copy()
                for ir in xrange(self.N_parproj[ishell]):
                    for sig,gf in tmp: tmp[sig] <<= self.downfold_pc(ik,ir,ishell,sig,S[sig],gf)
                    Gproj += tmp
                   
                # TO BE FIXED:
                # rotate to local frame
                #if (self.use_rotations):
                #    for sig,gf in Gproj: Gproj[sig] <<= self.rotloc(0,gf,direction='toLocal')

                for iom in range(N_om): 
                    if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                        for ish in range(self.shells[ishell][3]):
                            for ibn in bln:
                                Akw[ibn][ish,ik,iom] = Gproj[ibn]._data.array[ish,ish,iom].imag/(-3.1415926535)
               
            
        # END k-LOOP
        if (MPI.IS_MASTER_NODE()):
            if (ishell is None):
        
                for ibn in bln:
                    # loop over GF blocs:
    
                    if (invertAkw):
                        maxAkw=Akw[ibn].max()
                        minAkw=Akw[ibn].min()


                    # open file for storage:
                    if Fermisurface:
                        f=open('FS_'+ibn+'.dat','w')
                    else:
                        f=open('Akw_'+ibn+'.dat','w')

                    for ik in range(self.Nk):
                        if Fermisurface:
                            if (invertAkw):
                                Akw[ibn][ik,0] = 1.0/(minAkw-maxAkw)*(Akw[ibn][ik,0] - maxAkw)           
                            f.write('%s    %s\n'%(ik,Akw[ibn][ik,0]))
                        else:
                            for iom in range(N_om): 
                                if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                                    if (invertAkw):
                                        Akw[ibn][ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ibn][ik,iom] - maxAkw)
                                    if (shift>0.0001):
                                        f.write('%s      %s\n'%(M[iom],Akw[ibn][ik,iom]))
                                    else:
                                        f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ibn][ik,iom]))

                            f.write('\n')
 
                    f.close()

            else:
                for ibn in bln:
                    for ish in range(self.shells[ishell][3]):
                     
                        if (invertAkw):
                            maxAkw=Akw[ibn][ish,:,:].max()
                            minAkw=Akw[ibn][ish,:,:].min()

                        f=open('Akw_'+ibn+'_proj'+str(ish)+'.dat','w') 

                        for ik in range(self.Nk):
                            for iom in range(N_om): 
                                if (M[iom]>omminplot) and (M[iom]<ommaxplot):
                                    if (invertAkw):
                                        Akw[ibn][ish,ik,iom] = 1.0/(minAkw-maxAkw)*(Akw[ibn][ish,ik,iom] - maxAkw)
                                    if (shift>0.0001):
                                        f.write('%s      %s\n'%(M[iom],Akw[ibn][ish,ik,iom]))
                                    else:
                                        f.write('%s     %s      %s\n'%(ik,M[iom],Akw[ibn][ish,ik,iom]))

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


        

    def partial_charges(self):
        """Calculates the orbitally-resolved density matrix for all the orbitals considered in the input.
           The theta-projectors are used, hence case.parproj data is necessary"""
           

        #thingstoread = ['Dens_Mat_below','N_parproj','Proj_Mat_pc','rotmat_all']
        #retval = self.read_input_from_HDF(SubGrp=self.ParProjdata,thingstoread=thingstoread)
        retval = self.read_ParProj_input_from_HDF()
        if not retval: return retval
        if self.symm_op: self.Symm_par = Symmetry(self.HDFfile,subgroup=self.Symmpardata)
        
        # Density matrix in the window
        bln = self.blocnames[self.SO]
        ntoi = self.names_to_ind[self.SO]
        self.Dens_Mat_window = [ [numpy.zeros([self.shells[ish][3],self.shells[ish][3]],numpy.complex_) for ish in range(self.N_shells)]   
                                 for isp in range(len(bln)) ]    # init the density matrix

        mu = self.Chemical_Potential
        GFStruct_proj = [ [ (al, range(self.shells[i][3])) for al in bln ]  for i in xrange(self.N_shells) ]
        if hasattr(self,"Sigmaimp"):
            Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Mesh = self.Sigmaimp[0].mesh)) for a,al in GFStruct_proj[ish] ], Copy = False)
                     for ish in xrange(self.N_shells)]
        else:
            Gproj = [GF(Name_Block_Generator = [ (a,GFBloc_ImFreq(Indices = al, Beta = 40)) for a,al in GFStruct_proj[ish] ], Copy = False)
                     for ish in xrange(self.N_shells)]

        for ish in xrange(self.N_shells): Gproj[ish].zero()

        ikarray=numpy.array(range(self.Nk))
        #print MPI.rank, MPI.slice_array(ikarray)
        #print "K-Sum starts on node",MPI.rank," at ",datetime.now()
        
        for ik in MPI.slice_array(ikarray):
            #print MPI.rank, ik, datetime.now()
            S = self.latticeGF_Matsubara(ik=ik,mu=mu)
            S *= self.BZ_weights[ik]

            for ish in xrange(self.N_shells):
                tmp = Gproj[ish].copy()
                for ir in xrange(self.N_parproj[ish]):
                    for sig,gf in tmp: tmp[sig] <<= self.downfold_pc(ik,ir,ish,sig,S[sig],gf)
                    Gproj[ish] += tmp
        
        #print "K-Sum done on node",MPI.rank," at ",datetime.now()
        #collect data from MPI:
        for ish in xrange(self.N_shells):
            Gproj[ish] <<= MPI.all_reduce(MPI.world,Gproj[ish],lambda x,y : x+y)
        MPI.barrier()

        #print "Data collected on node",MPI.rank," at ",datetime.now()

        # Symmetrisation:
        if (self.symm_op!=0): Gproj = self.Symm_par.symmetrise(Gproj)
        #print "Symmetrisation done on node",MPI.rank," at ",datetime.now()
        
        for ish in xrange(self.N_shells):

            # Rotation to local:
            if (self.use_rotations):
                for sig,gf in Gproj[ish]: Gproj[ish][sig] <<= self.rotloc_all(ish,gf,direction='toLocal')

            isp = 0
            for sig,gf in Gproj[ish]: #dmg.append(Gproj[ish].density()[sig])
                self.Dens_Mat_window[isp][ish] = Gproj[ish].density()[sig]
                isp+=1
       
        # add Density matrices to get the total:
        Dens_Mat = [ [ self.Dens_Mat_below[ntoi[bln[isp]]][ish]+self.Dens_Mat_window[isp][ish] for ish in range(self.N_shells)]
                     for isp in range(len(bln)) ]

        return Dens_Mat



