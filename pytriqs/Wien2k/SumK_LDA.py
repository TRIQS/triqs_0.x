
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
from pytriqs.Wien2k.Symmetry import *
import numpy
import pytriqs.Base.Utility.Dichotomy as Dichotomy
from pytriqs.Base.GF_Local.GF import GF
from pytriqs.Base.GF_Local.GFBloc_ImFreq import GFBloc_ImFreq
from pytriqs.Base.GF_Local.GFBloc_ReFreq import GFBloc_ReFreq
from pytriqs.Base.GF_Local import GF_Initializers
from pytriqs.Solvers.Operators import *
from pytriqs.Base.Archive.HDF_Archive import *
import pytriqs.Base.Utility.MPI as MPI

from math import cos,sin

import string, pickle

class SumK_LDA:
    """This class provides a general SumK method for combining ab-initio code and pytriqs."""


    def __init__(self, HDFfile, mu = 0.0, hfield = 0.0, UseLDABlocs = False, LDAdata = 'SumK_LDA', Symmcorrdata = 'SymmCorr',
                 ParProjdata = 'SumK_LDA_ParProj', Symmpardata = 'SymmPar', Bandsdata = 'SumK_LDA_Bands'):
        """
        Initialises the class from data previously stored into an HDF5
        """

        if  not (type(HDFfile)==StringType):
            MPI.report("Give a string for the HDF5 filename to read the input!")
        else:
            self.HDFfile = HDFfile
            self.LDAdata = LDAdata
            self.ParProjdata = ParProjdata
            self.Bandsdata = Bandsdata
            self.Symmpardata = Symmpardata
            self.Symmcorrdata = Symmcorrdata
            self.blocnames= [ ['up','down'], ['ud'] ]
            self.NspinblocsGF = [2,1]
            self.Gupf = None
            self.hfield = hfield
            
            # read input from HDF:
            thingstoread = ['EnergyUnit','Nk','k_dep_projection','SP','SO','charge_below','Density_Required',
                            'symm_op','N_shells','shells','N_corr_shells','corr_shells','use_rotations','rotmat','rotmat_timeinv','Nreps',
                            'dim_reps','T','N_Orbitals','Proj_Mat','BZ_weights','Hopping']
            optionalthings = ['GFStruct_Solver','mapinv','map','Chemical_Potential','dc_imp','DCenerg','deg_shells']

            #ar=HDF_Archive(self.HDFfile,'a')
            #del ar

            self.retval = self.read_input_from_HDF(SubGrp=self.LDAdata,thingstoread=thingstoread,optionalthings=optionalthings)

            #ar=HDF_Archive(self.HDFfile,'a')
            #del ar

            if (self.SO) and (abs(self.hfield)>0.000001):
                self.hfield=0.0
                MPI.report("For SO, the external magnetic field is not implemented, setting it to 0!!")

           
            self.inequiv_shells(self.corr_shells)     # determine the number of inequivalent correlated shells

            # field to convert blocnames to indices
            self.names_to_ind = [{}, {}]
            for ibl in range(2):
                for inm in range(self.NspinblocsGF[ibl]): 
                    self.names_to_ind[ibl][self.blocnames[ibl][inm]] = inm * self.SP #(self.Nspinblocs-1)

            # GF structure used for the local things in the k sums
            self.GFStruct_corr = [ [ (al, range( self.corr_shells[i][3])) for al in self.blocnames[self.corr_shells[i][4]] ]  
                                   for i in xrange(self.N_corr_shells) ]

            if not (self.retval['GFStruct_Solver']):
                # No GFStruct was stored in HDF, so first set a standard one:
                self.GFStruct_Solver = [ [ (al, range( self.corr_shells[self.invshellmap[i]][3]) )
                                           for al in self.blocnames[self.corr_shells[self.invshellmap[i]][4]] ]
                                         for i in xrange(self.N_inequiv_corr_shells) ]
                self.map = [ {} for i in xrange(self.N_inequiv_corr_shells) ]
                self.mapinv = [ {} for i in xrange(self.N_inequiv_corr_shells) ]
                for i in xrange(self.N_inequiv_corr_shells):
                    for al in self.blocnames[self.corr_shells[self.invshellmap[i]][4]]:
                        self.map[i][al] = [al for j in range( self.corr_shells[self.invshellmap[i]][3] ) ]
                        self.mapinv[i][al] = al

            if not (self.retval['dc_imp']):
                # init the double counting:
                self.__initDC()

            if not (self.retval['Chemical_Potential']):
                self.Chemical_Potential = mu

            if not (self.retval['deg_shells']):
                self.deg_shells = [ [] for i in range(self.N_inequiv_corr_shells)]

            if self.symm_op:
                #MPI.report("Do the init for symm:")
                self.Symm_corr = Symmetry(HDFfile,subgroup=self.Symmcorrdata)

            # determine the smallest blocs, if wanted:
            if (UseLDABlocs): dm=self.analyse_BS()

          
            # now save things again to HDF5:
            if (MPI.IS_MASTER_NODE()):
                ar=HDF_Archive(self.HDFfile,'a')
                ar[self.LDAdata]['hfield'] = self.hfield
                del ar
            self.save()

            

            
 

            
    def read_input_from_HDF(self,SubGrp,thingstoread,optionalthings=[]):
        """
        Reads data from the HDF file
        """
        
        retval = True
        # init variables on all nodes:
        for it in thingstoread: exec "self.%s = 0"%it
        for it in optionalthings: exec "self.%s = 0"%it
        
        if (MPI.IS_MASTER_NODE()):
            ar=HDF_Archive(self.HDFfile,'a')
            if (SubGrp in ar):
                # first read the necessary things:
                for it in thingstoread:
                    if (it in ar[SubGrp]):
                        exec "self.%s = ar['%s']['%s']"%(it,SubGrp,it)
                    else:
                        MPI.report("Loading %s failed!"%it)
                        retval = False
                   
                if ((retval) and (len(optionalthings)>0)):
                    # if necessary things worked, now read optional things:
                    retval = {}
                    for it in optionalthings:
                        if (it in ar[SubGrp]):
                            exec "self.%s = ar['%s']['%s']"%(it,SubGrp,it)
                            retval['%s'%it] = True
                        else:
                            retval['%s'%it] = False
            else:
                MPI.report("Loading failed: No %s subgroup in HDF5!"%SubGrp)
                retval = False

            del ar

        # now do the broadcasting:
        for it in thingstoread: exec "self.%s = MPI.bcast(self.%s)"%(it,it)
        for it in optionalthings: exec "self.%s = MPI.bcast(self.%s)"%(it,it)
        

        retval = MPI.bcast(retval)
               
        return retval



    def save(self):
        """Saves some quantities into an HDF5 arxiv"""

        if not (MPI.IS_MASTER_NODE()): return # do nothing on nodes

        ar=HDF_Archive(self.HDFfile,'a')
        ar[self.LDAdata]['Chemical_Potential'] = self.Chemical_Potential
        ar[self.LDAdata]['DCenerg'] = self.DCenerg
        ar[self.LDAdata]['dc_imp'] = self.dc_imp
        del ar      
            

    def load(self):
        """Loads some quantities from an HDF5 arxiv"""

        thingstoread=['Chemical_Potential','dc_imp','DCenerg']
        
        retval = self.read_input_from_HDF(SubGrp=self.LDAdata,thingstoread=thingstoread)
        return retval


    def downfold(self,ik,icrsh,sig,gf_to_downfold,gf_inp):
        """Downfolding a block of the Greens function"""
        
        gf_downfolded = gf_inp.copy()
        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        gf_downfolded.from_L_G_R(self.Proj_Mat[ik][isp][icrsh],gf_to_downfold,self.Proj_Mat[ik][isp][icrsh].conjugate().transpose())  # downfolding G

        return gf_downfolded
        

    def upfold(self,ik,icrsh,sig,gf_to_upfold,gf_inp):
        """Upfolding a block of the Greens function"""

        gf_upfolded = gf_inp.copy()
        
        isp = self.names_to_ind[self.SO][sig]       # get spin index for proj. matrices
        gf_upfolded.from_L_G_R(self.Proj_Mat[ik][isp][icrsh].conjugate().transpose(),gf_to_upfold,self.Proj_Mat[ik][isp][icrsh]) 

        return gf_upfolded


    def rotloc(self,icrsh,gf_to_rotate,direction):
        """Local <-> Global rotation of a GF block.
           direction: 'toLocal' / 'toGlobal' """

        assert ((direction=='toLocal')or(direction=='toGlobal')),"Give direction 'toLocal' or 'toGlobal' in rotloc!"

        gf_rotated = gf_to_rotate.copy()
        if (direction=='toGlobal'):
            #if (self.rotmat_timeinv[icrsh]==1): gf_rotated <<= gf_rotated.transpose()
            #gf_rotated.from_L_G_R(self.rotmat[icrsh].transpose(),gf_rotated,self.rotmat[icrsh].conjugate())
            if ((self.rotmat_timeinv[icrsh]==1) and (self.SO)):
                gf_rotated <<= gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rotmat[icrsh].conjugate(),gf_rotated,self.rotmat[icrsh].transpose())
            else:
                gf_rotated.from_L_G_R(self.rotmat[icrsh],gf_rotated,self.rotmat[icrsh].conjugate().transpose())
            
        elif (direction=='toLocal'):
            if ((self.rotmat_timeinv[icrsh]==1)and(self.SO)):
                gf_rotated <<= gf_rotated.transpose()
                gf_rotated.from_L_G_R(self.rotmat[icrsh].transpose(),gf_rotated,self.rotmat[icrsh].conjugate())
            else:
                gf_rotated.from_L_G_R(self.rotmat[icrsh].conjugate().transpose(),gf_rotated,self.rotmat[icrsh])
                
        return gf_rotated
                    

    def latticeGF_Matsubara(self,ik,mu,Beta=40,withSigma=True):
        """Calculates the lattice Green function from the LDA hopping and the self energy at k-point number ik
           and chemical potential mu."""

        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        if (not hasattr(self,"Sigmaimp")): withSigma=False

        if (withSigma): 
            stmp = self.add_DC()
            Beta = self.Sigmaimp[0].Beta        #override Beta if Sigma is present
            
        if (self.Gupf is None):
            # first setting up of Gupf
            BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
            GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
            a_list = [a for a,al in GFStruct]   
            glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta = Beta) for a,al in GFStruct]  
            self.Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
            self.Gupf.zero()

        GFsize = [ gf.N1 for sig,gf in self.Gupf]  
        unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==GFsize[ib] 
                               for ib in range(self.NspinblocsGF[self.SO]) ] )

        if (not unchangedsize):
            BS = [ range(self.N_Orbitals[ik][ntoi[ib]]) for ib in bln ]
            GFStruct = [ (bln[ib], BS[ib]) for ib in range(self.NspinblocsGF[self.SO]) ]
            a_list = [a for a,al in GFStruct]                                 
            glist = lambda : [ GFBloc_ImFreq(Indices = al, Beta = Beta) for a,al in GFStruct]    
            self.Gupf = GF(NameList = a_list, BlockList = glist(),Copy=False)
            self.Gupf.zero()

        idmat = [numpy.identity(self.N_Orbitals[ik][ntoi[bl]],numpy.complex_) for bl in bln]  
        #for ibl in range(self.NspinblocsGF[self.SO]): mupat[ibl] *= mu

        self.Gupf <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
        M = copy.deepcopy(idmat)
        for ibl in range(self.NspinblocsGF[self.SO]): 
            ind = ntoi[bln[ibl]]
            M[ibl] = self.Hopping[ik][ind] - (idmat[ibl]*mu) - (idmat[ibl] * self.hfield * (1-2*ibl))
        self.Gupf -= M

        if (withSigma):
            tmp = self.Gupf.copy()    # init temporary storage
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in tmp: tmp[sig] <<= self.upfold(ik,icrsh,sig,stmp[icrsh][sig],gf)
                self.Gupf -= tmp      # adding to the upfolded GF

        self.Gupf.invert()

	return self.Gupf


    def check_projectors(self):

        densmat = [numpy.zeros([self.corr_shells[ish][3],self.corr_shells[ish][3]],numpy.complex_) 
                   for ish in range(self.N_corr_shells)]
        
        for ik in range(self.Nk):
        
            for ish in range(self.N_corr_shells):
                Norb = self.corr_shells[ish][3]
                densmat[ish][:,:] += numpy.dot(self.Proj_Mat[ik][0][ish],self.Proj_Mat[ik][0][ish].transpose().conjugate()) * self.BZ_weights[ik]

        if (self.symm_op!=0): densmat = self.Symm_corr.symmetrise(densmat)

        # Rotate to local coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                if (self.rotmat_timeinv[icrsh]==1): densmat[icrsh] = densmat[icrsh].conjugate()
                densmat[icrsh] = numpy.dot( numpy.dot(self.rotmat[icrsh].conjugate().transpose(),densmat[icrsh]) , 
                                            self.rotmat[icrsh] )
                
               
        return densmat



    def simplepointdensmat(self):


        ntoi = self.names_to_ind[self.SO]
        bln = self.blocnames[self.SO]

        MMat = [numpy.zeros( [self.N_Orbitals[0][ntoi[bl]],self.N_Orbitals[0][ntoi[bl]]], numpy.complex_) for bl in bln] 

        densmat = [ {} for icrsh in xrange(self.N_corr_shells)]
        for icrsh in xrange(self.N_corr_shells):
            for bl in self.blocnames[self.corr_shells[icrsh][4]]:
                densmat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]], numpy.complex_)

        ikarray=numpy.array(range(self.Nk))
          
        for ik in MPI.slice_array(ikarray):
            
            unchangedsize = all( [ self.N_Orbitals[ik][ntoi[bln[ib]]]==len(MMat[ib]) 
                                   for ib in range(self.NspinblocsGF[self.SO]) ] )
               
            if (not unchangedsize):
                MMat = [numpy.zeros( [self.N_Orbitals[ik][ntoi[bl]],self.N_Orbitals[ik][ntoi[bl]]], numpy.complex_) for bl in bln] 

            for ibl,bl in enumerate(bln):
                ind = ntoi[bl]
                for inu in range(self.N_Orbitals[ik][ind]):
                    if ( (self.Hopping[ik][ind][inu,inu]-self.hfield*(1-2*ibl)) < 0.0): 
                        MMat[ibl][inu,inu] = 1.0
                    else:
                        MMat[ibl][inu,inu] = 0.0 


            for icrsh in range(self.N_corr_shells):
                for ibn,bn in enumerate(self.blocnames[self.corr_shells[icrsh][4]]):
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    #print ik, bn, isp
                    densmat[icrsh][bn] += self.BZ_weights[ik] * numpy.dot( numpy.dot(self.Proj_Mat[ik][isp][icrsh],MMat[ibn]) , 
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
                    if (self.rotmat_timeinv[icrsh]==1): densmat[icrsh][bn] = densmat[icrsh][bn].conjugate()
                    densmat[icrsh][bn] = numpy.dot( numpy.dot(self.rotmat[icrsh].conjugate().transpose(),densmat[icrsh][bn]) , 
                                                    self.rotmat[icrsh])
                

        return densmat


    def density_gf(self,Beta = 40):
        """Calculates the density without setting up Gloc. It is useful for Hubbard I, and very fast.""" 

        densmat = [ {} for icrsh in xrange(self.N_corr_shells)]
        for icrsh in xrange(self.N_corr_shells):
            for bl in self.blocnames[self.corr_shells[icrsh][4]]:
                densmat[icrsh][bl] = numpy.zeros([self.corr_shells[icrsh][3],self.corr_shells[icrsh][3]], numpy.complex_)

        ikarray=numpy.array(range(self.Nk))

        for ik in MPI.slice_array(ikarray):
            
            Gupf = self.latticeGF_Matsubara(ik=ik,mu=self.Chemical_Potential)
            Gupf *= self.BZ_weights[ik]
            dm = Gupf.density()
            MMat = [dm[bl] for bl in self.blocnames[self.SO]]
            
            for icrsh in range(self.N_corr_shells):
                for ibn,bn in enumerate(self.blocnames[self.corr_shells[icrsh][4]]):
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    #print ik, bn, isp
                    densmat[icrsh][bn] += numpy.dot( numpy.dot(self.Proj_Mat[ik][isp][icrsh],MMat[ibn]),self.Proj_Mat[ik][isp][icrsh].transpose().conjugate() )

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
                    if (self.rotmat_timeinv[icrsh]==1): densmat[icrsh][bn] = densmat[icrsh][bn].conjugate()
                    densmat[icrsh][bn] = numpy.dot( numpy.dot(self.rotmat[icrsh].conjugate().transpose(),densmat[icrsh][bn]) , 
                                                    self.rotmat[icrsh] )

        return densmat



    def analyse_BS(self, threshold = 0.00001, includeshells = None):
        """ Determines the Greens function block structure from the simple point integration"""

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


            # now calculate degeneracies of orbitals:
            dm = {}
            for bl in self.GFStruct_Solver[ish]:
                bln = bl[0]
                ind = bl[1]
                # get dm for the blocks:
                dm[bln] = numpy.zeros([len(ind),len(ind)],numpy.complex_)
                for i in range(len(ind)):
                    for j in range(len(ind)):
                        dm[bln][i,j] = densmat[ish][self.mapinv[ish][bln]][ind[i],ind[j]]

            for bl in self.GFStruct_Solver[ish]:
                for bl2 in self.GFStruct_Solver[ish]:
                    if (dm[bl[0]].shape==dm[bl2[0]].shape) :
                        if ( ( (abs(dm[bl[0]]-dm[bl2[0]])<threshold).all() ) and (bl[0]!=bl2[0]) ):
                            # check if it was already there:
                            ind1=-1
                            ind2=-2
                            for n,ind in enumerate(self.deg_shells[ish]):
                                if (bl[0] in ind): ind1=n
                                if (bl2[0] in ind): ind2=n
                            if ((ind1<0)and(ind2>=0)):
                                self.deg_shells[ish][ind2].append(bl[0])
                            elif ((ind1>=0)and(ind2<0)):
                                self.deg_shells[ish][ind1].append(bl2[0])
                            elif ((ind1<0)and(ind2<0)):
                                self.deg_shells[ish].append([bl[0],bl2[0]])

        if (MPI.IS_MASTER_NODE()):
            ar=HDF_Archive(self.HDFfile,'a')
            ar[self.LDAdata]['GFStruct_Solver'] = self.GFStruct_Solver
            ar[self.LDAdata]['map'] = self.map
            ar[self.LDAdata]['mapinv'] = self.mapinv
            try:
                ar[self.LDAdata]['deg_shells'] = self.deg_shells
            except:
                MPI.report("deg_shells not stored, degeneracies not found")
            del ar
            
        return densmat
        

    def symm_deg_GF(self,gftosymm,orb):
        """Symmetrises a GF for the given degenerate shells self.deg_shells"""
        
        for degsh in self.deg_shells[orb]:
            #loop over degenerate shells:
            ss = gftosymm[degsh[0]].copy()
            ss.zero()
            Ndeg = len(degsh)
            for bl in degsh: ss += gftosymm[bl] / (1.0*Ndeg)
            for bl in degsh: gftosymm[bl] <<= ss
        

    def eff_atomic_levels(self):
        """Calculates the effective atomic levels needed as input for the Hubbard I Solver."""

        # define matrices for inequivalent shells:
        eff_atlevels = [ {} for ish in range(self.N_inequiv_corr_shells) ]
        for ish in range(self.N_inequiv_corr_shells):
            for bn in self.blocnames[self.corr_shells[self.invshellmap[ish]][4]]:
                eff_atlevels[ish][bn] = numpy.identity(self.corr_shells[self.invshellmap[ish]][3], numpy.complex_)
 
        # Chemical Potential:
        for ish in xrange(self.N_inequiv_corr_shells): 
            for ii in eff_atlevels[ish]: eff_atlevels[ish][ii] *= -self.Chemical_Potential
        
        # double counting term:
        #if hasattr(self,"dc_imp"):
        for ish in xrange(self.N_inequiv_corr_shells): 
            for ii in eff_atlevels[ish]:
                eff_atlevels[ish][ii] -= self.dc_imp[self.invshellmap[ish]][ii]

        # sum over k:
        if not hasattr(self,"Hsumk"):
            # calculate the sum over k. Does not depend on mu, so do it only once:
            self.Hsumk = [ {} for ish in range(self.N_corr_shells) ]
            for icrsh in range(self.N_corr_shells):
                for bn in self.blocnames[self.corr_shells[icrsh][4]]: 
                    dim = self.corr_shells[icrsh][3]  #*(1+self.corr_shells[icrsh][4])
                    self.Hsumk[icrsh][bn] = numpy.zeros([dim,dim],numpy.complex_)
          
            for icrsh in range(self.N_corr_shells):
                for ibn, bn in enumerate(self.blocnames[self.corr_shells[icrsh][4]]):
                    isp = self.names_to_ind[self.corr_shells[icrsh][4]][bn]
                    for ik in xrange(self.Nk):
                        MMat = numpy.identity(self.N_Orbitals[ik][isp], numpy.complex_)
                        MMat = self.Hopping[ik][isp] - (1-2*ibn) * self.hfield * MMat
                        self.Hsumk[icrsh][bn] += self.BZ_weights[ik] * numpy.dot( numpy.dot(self.Proj_Mat[ik][isp][icrsh],MMat), #self.Hopping[ik][isp]) , 
                                                                                  self.Proj_Mat[ik][isp][icrsh].conjugate().transpose() )

            # symmetrisation:
            if (self.symm_op!=0): self.Hsumk = self.Symm_corr.symmetrise(self.Hsumk)

            # Rotate to local coordinate system:
            if (self.use_rotations):
                for icrsh in xrange(self.N_corr_shells):
                    for bn in self.Hsumk[icrsh]:

                        if (self.corr_shells[icrsh][4]==0): self.Hsumk[icrsh][bn] = self.Hsumk[icrsh][bn].conjugate()
                        self.Hsumk[icrsh][bn] = numpy.dot( numpy.dot(self.rotmat[icrsh].conjugate().transpose(),self.Hsumk[icrsh][bn]) , 
                                                           self.rotmat[icrsh] )
                 
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
            l = self.corr_shells[i][3]
            for j in xrange(len(self.GFStruct_corr[i])):
                self.dc_imp[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.zeros([l,l],numpy.float_)
        self.DCenerg = [0.0 for i in xrange(self.N_corr_shells)]



    def SetDC_Lichtenstein(self,sigimp):
        """Sets a double counting term according to Lichtenstein et al. PRL2001"""

        assert isinstance(sigimp,list), "sigimp has to be a list of impurity self energies for the correlated shells, even if it is of length 1!"
        assert len(sigimp)==self.N_inequiv_corr_shells, "give exactly one Sigma for each inequivalent corr. shell!"

        for i in xrange(self.N_corr_shells):
            l = (self.corr_shells[i][4]+1) * self.corr_shells[i][3]
            for j in xrange(len(self.GFStruct_corr[i])):
                self.dc_imp[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.identity(l,numpy.float_)
 
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
        

        #if (not hasattr(self,"dc_imp")): self.__initDC()
                    
                
        dm = [ {} for i in xrange(self.N_corr_shells)]
        for i in xrange(self.N_corr_shells):
            l = self.corr_shells[i][3] #*(1+self.corr_shells[i][4])
            for j in xrange(len(self.GFStruct_corr[i])):
                dm[i]['%s'%self.GFStruct_corr[i][j][0]] = numpy.zeros([l,l],numpy.float_)
        

        for icrsh in xrange(self.N_corr_shells):

            iorb = self.shellmap[icrsh]    # iorb is the index of the inequivalent shell corresponding to icrsh

            if (iorb==orb):
                # do this orbital

                l = self.corr_shells[icrsh][3] #*(1+self.corr_shells[icrsh][4])
                for j in xrange(len(self.GFStruct_corr[icrsh])):
                    self.dc_imp[icrsh]['%s'%self.GFStruct_corr[icrsh][j][0]] = numpy.identity(l,numpy.float_)


                # transform the CTQMC blocks to the full matrix:
                for ibl in range(len(self.GFStruct_Solver[iorb])):
                    for i in range(len(self.GFStruct_Solver[iorb][ibl][1])):
                        for j in range(len(self.GFStruct_Solver[iorb][ibl][1])):
                            bl   = self.GFStruct_Solver[iorb][ibl][0]
                            ind1 = self.GFStruct_Solver[iorb][ibl][1][i]
                            ind2 = self.GFStruct_Solver[iorb][ibl][1][j]
                            dm[icrsh][self.mapinv[iorb][bl]][ind1,ind2] = densmat[bl][i,j].real    # only real part relevant for trace

                M = self.corr_shells[icrsh][3]
                Ncr = {}
                Ncrtot = 0.0
                a_list = [a for a,al in self.GFStruct_corr[icrsh]]
                for bl in a_list:
                    Ncr[bl] = dm[icrsh][bl].trace()
                    Ncrtot += Ncr[bl]

                # average the densities if there is no SP:
                if (self.SP==0):
                    for bl in a_list:
                        Ncr[bl] = Ncrtot / len(a_list)
                

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
        """Puts the impurity self energies for inequivalent atoms into the class, respects the multiplicity of the atoms."""

        assert isinstance(Sigmaimp,list), "Sigmaimp has to be a list of Sigmas for the correlated shells, even if it is of length 1!"
        assert len(Sigmaimp)==self.N_inequiv_corr_shells, "give exactly one Sigma for each inequivalent corr. shell!"

       
        # init self.Sigmaimp:
        if (Sigmaimp[0].Note=='ReFreq'):
            # Real frequency Sigma:
            self.Sigmaimp = [ GF( Name_Block_Generator = [ (a,GFBloc_ReFreq(Indices = al, Mesh = Sigmaimp[0].mesh)) for a,al in self.GFStruct_corr[i] ],
                                  Copy = False) for i in xrange(self.N_corr_shells) ]
            self.Sigmaimp[0].Note = 'ReFreq'
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

        # rotation from local to global coordinate system:
        if (self.use_rotations):
            for icrsh in xrange(self.N_corr_shells):
                for sig,gf in self.Sigmaimp[icrsh]: self.Sigmaimp[icrsh][sig] <<= self.rotloc(icrsh,gf,direction='toGlobal')



    def add_DC(self):
        """Substracts the double counting term from the impurity self energy."""
        
        # Be careful: Sigmaimp is already in the global coordinate system!!
        sres = [s.copy() for s in self.Sigmaimp]
        for icrsh in xrange(self.N_corr_shells):
            for bl,gf in sres[icrsh]:
                dccont = numpy.dot(self.rotmat[icrsh],numpy.dot(self.dc_imp[icrsh][bl],self.rotmat[icrsh].conjugate().transpose()))
                sres[icrsh][bl] -= dccont
            
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
        
        The calculation is done in the global coordinate system, if distinction is made between local/global!
        """
        
        dens = 0.0
        ikarray=numpy.array(range(self.Nk))
        
        for ik in MPI.slice_array(ikarray):
        
            S = self.latticeGF_Matsubara(ik=ik,mu=mu) 
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
            
        Gloc = [ self.Sigmaimp[icrsh].copy() for icrsh in xrange(self.N_corr_shells) ]   # this list will be returned  
        for icrsh in xrange(self.N_corr_shells): Gloc[icrsh].zero()                # initialize to zero

        ikarray=numpy.array(range(self.Nk))
        
        for ik in MPI.slice_array(ikarray):
            
            S = self.latticeGF_Matsubara(ik=ik,mu=mu,withSigma = withSigma) 
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

        ikarray=numpy.array(range(self.Nk))
        
        dens = 0.0
        
        for ik in MPI.slice_array(ikarray):
        
            S = self.latticeGF_Matsubara(ik=ik,mu=self.Chemical_Potential)
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
                            valre = (deltaN['up'][ik][inu,imu].real + deltaN['down'][ik][inu,imu].real) / 2.0
                            valim = (deltaN['up'][ik][inu,imu].imag + deltaN['down'][ik][inu,imu].imag) / 2.0
                            f.write("%.14f  %.14f "%(valre,valim))
                        f.write("\n")
                    f.write("\n")
                f.close()
            elif ((self.SP==1)and(self.SO==0)):
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
            else:
                for ik in range(self.Nk):
                    f.write("%s\n"%self.N_Orbitals[ik][0])
                    for inu in range(self.N_Orbitals[ik][0]):
                        for imu in range(self.N_Orbitals[ik][0]):
                            f.write("%.14f  %.14f "%(deltaN['ud'][ik][inu,imu].real,deltaN['ud'][ik][inu,imu].imag))
                        f.write("\n")
                    f.write("\n")
                f.close()
                for ik in range(self.Nk):
                    f1.write("%s\n"%self.N_Orbitals[ik][0])
                    for inu in range(self.N_Orbitals[ik][0]):
                        for imu in range(self.N_Orbitals[ik][0]):
                            f1.write("%.14f  %.14f "%(deltaN['ud'][ik][inu,imu].real,deltaN['ud'][ik][inu,imu].imag))
                        f1.write("\n")
                    f1.write("\n")
                f1.close()
                                                            

        return deltaN, dens
