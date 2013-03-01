
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


from pytriqs.dft.U_matrix import *
from pytriqs.solvers.operators import *
from pytriqs.solvers.ctqmc_hyb import Solver
from pytriqs.base.utility.my_utils import sum_list
import pytriqs.base.utility.mpi as mpi
from types import *
import numpy


#########################################
#
#  Solver for the Multi-Band problem
#
#########################################


class SolverMultiBand (Solver):
    """ 
    This is a general solver for a multiband local Hamiltonian. 
    Calling arguments: 
    beta = inverse temperature
    n_orb = Number of local orbitals
    U_interact = Average Coulomb interaction
    J_hund     = Hund coupling
    use_spinflip = true/false
    use_pairhop  = true/false
    use_matrix: Use the interaction matrix calculated from the Slater integrals
    is use_matrix, you need also:
        l: angular momentum of the orbital, l=2 is d
        T: Transformation matrix for U vertex. If not present, use standard complex harmonics
               
    """
    
    def __init__(self, beta, n_orb, U_interact=None, J_hund=None, gf_struct=False, map=False, use_spinflip=False,
                 use_matrix = True, l=2, T=None, dim_reps=None, irep=None, deg_orbs = [], sl_int = None):
    
        self.offset = 0
        self.use_spinflip = use_spinflip
        self.n_orb = n_orb
       
        self.U, self.Up, self.U4ind, self.offset = set_U_matrix(U_interact,J_hund,n_orb,l,use_matrix,T,sl_int,use_spinflip,dim_reps,irep) 

        if (gf_struct):
            assert map, "give also the mapping!"
            self.map = map
        else:
            # standard gf_struct and map
            gf_struct = [ ('%s'%(ud),[n for n in range(n_orb)]) for ud in ['up','down'] ]
            self.map = {'up' : ['up' for v in range(self.n_orb)], 'down' : ['down' for v in range(self.n_orb)]}

        #print gf_struct,self.map
        
        if (use_spinflip==False):
            Hamiltonian = self.__set_hamiltonian_density()
        else:
            if (use_matrix):
                #Hamiltonian = self.__set_full_hamiltonian_slater()
                Hamiltonian = self.__set_spinflip_hamiltonian_slater()
            else:
                Hamiltonian = self.__set_full_hamiltonian_kanamori(J_hund = J_hund)

        Quantum_Numbers = self.__set_quantum_numbers(gf_struct)
    
        # Determine if there are only blocs of size 1:
        self.blocssizeone = True
        for ib in gf_struct:
            if (len(ib[1])>1): self.blocssizeone = False

       
        # now initialize the solver with the stuff given above:
        Solver.__init__(self,
                        Beta = beta,
                        GFstruct = gf_struct,
                        H_Local = Hamiltonian,
                        Quantum_Numbers = Quantum_Numbers )

        #self.set_global_moves(deg_orbs)

        self.N_Cycles  = 10000
        self.Nmax_Matrix = 100
        self.N_Time_Slices_Delta= 10000
        #if ((len(gf_struct)==2*n_orb) and (use_spinflip==False)): 
        if ((self.blocssizeone) and (use_spinflip==False)):
            self.Use_Segment_Picture = True
        else:
            self.Use_Segment_Picture = False
        # check what all these parameters do!!!
   

    def set_global_moves(self, deg_orbs, factor=0.05):
        # Sets some global moves given orbital degeneracies:
        
        strbl  = ''
        strind = ''
        inddone = []

        for orbs in deg_orbs:
            ln = len(orbs)
            orbsorted = sorted(orbs)
            for ii in range(ln):
                if (strbl!=''): strbl += ','
                bl1 = orbsorted[ii]
                bl2 = orbsorted[(ii+1)%ln]
                ind1 = [ll for ll in self.Sigma[bl1].indices ]
                ind2 = [ll for ll in self.Sigma[bl2].indices ]

                strbl += "'" + bl1 + "':'" + bl2 + "'"
                for kk, ind in enumerate(ind1):
                    if not (ind in inddone):
                        if (strind!=''): strind += ','
                        strind += '%s:%s'%(ind1[kk],ind2[kk])
                        inddone.append(ind)
                

        if len(deg_orbs)>0:
            str = 'self.Global_Moves = [ (%s, lambda (a,alpha,dag) : ({ '%factor + strbl + ' }[a], {' + strind + '}[alpha], dag) )]'
            exec str
    
        

    def __set_hamiltonian_density(self):
        # density-density Hamiltonian:
        
        spinblocs = [v for v in self.map]
        #print spinblocs
        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it

        for sp1 in spinblocs:
            for sp2 in spinblocs:
                for i in range(self.n_orb):
                    for j in range(self.n_orb):
                        if (sp1==sp2):
                            Hamiltonian += 0.5 * self.U[self.offset+i,self.offset+j] * N(self.map[sp1][i],i) * N(self.map[sp2][j],j) 
                        else:
                            Hamiltonian += 0.5 * self.Up[self.offset+i,self.offset+j] * N(self.map[sp1][i],i) * N(self.map[sp2][j],j) 

        Hamiltonian -= N(self.map[spinblocs[0]][0],0)      # substract the initializing value

        return Hamiltonian


    def __set_full_hamiltonian_slater(self):
      
        spinblocs = [v for v in self.map]
        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it
        #print "Starting..."
        # use the full 4-index U-matrix:
        #for sp1 in spinblocs:
        #    for sp2 in spinblocs:
        for m1 in range(self.n_orb):
            for m2 in range(self.n_orb):
                for m3 in range(self.n_orb):
                    for m4 in range(self.n_orb):
                        if (abs(self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4])>0.00001):
                            for sp1 in spinblocs:
                                for sp2 in spinblocs:
                                    #print sp1,sp2,m1,m2,m3,m4
                                    Hamiltonian += 0.5 * self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4] * \
                                        Cdag(self.map[sp1][m1],m1) * Cdag(self.map[sp2][m2],m2) * C(self.map[sp2][m4],m4) * C(self.map[sp1][m3],m3)
        #print "end..."
        Hamiltonian -= N(self.map[spinblocs[0]][0],0)      # substract the initializing value
                        
        return Hamiltonian


    def __set_spinflip_hamiltonian_slater(self):
        """Takes only spin-flip and pair-hopping terms"""
        
        spinblocs = [v for v in self.map]
        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it
        #print "Starting..."
        # use the full 4-index U-matrix:
        #for sp1 in spinblocs:
        #    for sp2 in spinblocs:
        for m1 in range(self.n_orb):
            for m2 in range(self.n_orb):
                for m3 in range(self.n_orb):
                    for m4 in range(self.n_orb):
                        if ((abs(self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4])>0.00001) and
                            ( ((m1==m2)and(m3==m4)) or ((m1==m3)and(m2==m4)) or ((m1==m4)and(m2==m3)) ) ):
                            for sp1 in spinblocs:
                                for sp2 in spinblocs:
                                    #print sp1,sp2,m1,m2,m3,m4
                                    Hamiltonian += 0.5 * self.U4ind[self.offset+m1,self.offset+m2,self.offset+m3,self.offset+m4] * \
                                        Cdag(self.map[sp1][m1],m1) * Cdag(self.map[sp2][m2],m2) * C(self.map[sp2][m4],m4) * C(self.map[sp1][m3],m3)
        #print "end..."
        Hamiltonian -= N(self.map[spinblocs[0]][0],0)      # substract the initializing value
                        
        return Hamiltonian


            
    def __set_full_hamiltonian_kanamori(self,J_hund):

        spinblocs = [v for v in self.map]
        assert len(spinblocs)==2,"spinflips in Kanamori representation only implemented for up/down structure!"

        Hamiltonian = N(self.map[spinblocs[0]][0],0)       # initialize it

        # density terms:
        for sp1 in spinblocs:
            for sp2 in spinblocs:
                for i in range(self.n_orb):
                    for j in range(self.n_orb):
                        if (sp1==sp2):
                            Hamiltonian += 0.5 * self.U[self.offset+i,self.offset+j] * N(self.map[sp1][i],i) * N(self.map[sp2][j],j) 
                        else: 
                            Hamiltonian += 0.5 * self.Up[self.offset+i,self.offset+j] * N(self.map[sp1][i],i) * N(self.map[sp2][j],j) 

        # spinflip term:
        sp1 = spinblocs[0]
        sp2 = spinblocs[1]
        for i in range(self.n_orb-1):
            for j in range(i+1,self.n_orb):
                Hamiltonian -= J_hund * ( Cdag(self.map[sp1][i],i) * C(self.map[sp2][i],i) * Cdag(self.map[sp2][j],j) * C(self.map[sp1][j],j) )     # first term
                Hamiltonian -= J_hund * ( Cdag(self.map[sp2][i],i) * C(self.map[sp1][i],i) * Cdag(self.map[sp1][j],j) * C(self.map[sp2][j],j) )     # second term

        # pairhop terms:
        for i in range(self.n_orb-1):
            for j in range(i+1,self.n_orb):
                Hamiltonian -= J_hund * ( Cdag(self.map[sp1][i],i) * Cdag(self.map[sp2][i],i) * C(self.map[sp1][j],j) * C(self.map[sp2][j],j) )     # first term
                Hamiltonian -= J_hund * ( Cdag(self.map[sp2][j],j) * Cdag(self.map[sp1][j],j) * C(self.map[sp2][i],i) * C(self.map[sp1][i],i) )     # second term  

        Hamiltonian -= N(self.map[spinblocs[0]][0],0)       # substract the initializing value
                        
        return Hamiltonian
   

    def __set_quantum_numbers(self,gf_struct):
    
        QN = {}
        spinblocs = [v for v in self.map]

        # Define the quantum numbers:
        if (self.use_spinflip) :            
            Ntot = sum_list( [ N(self.map[s][i],i) for s in spinblocs for i in range(self.n_orb) ] )
            QN['NtotQN'] = Ntot
            #QN['Ntot'] = sum_list( [ N(self.map[s][i],i) for s in spinblocs for i in range(self.n_orb) ] )
            if (len(spinblocs)==2):
                # Assuming up/down structure:
                Sz = sum_list( [ N(self.map[spinblocs[0]][i],i)-N(self.map[spinblocs[1]][i],i) for i in range(self.n_orb) ] )
                QN['SzQN'] = Sz
                # new quantum number: works only if there are only spin-flip and pair hopping, not any more complicated things
                for i in range(self.n_orb):
                    QN['Sz2_%s'%i] = (N(self.map[spinblocs[0]][i],i)-N(self.map[spinblocs[1]][i],i)) * (N(self.map[spinblocs[0]][i],i)-N(self.map[spinblocs[1]][i],i))

        else :
            for ibl in range(len(gf_struct)):
                QN['N%s'%gf_struct[ibl][0]] = sum_list( [ N(gf_struct[ibl][0],gf_struct[ibl][1][i]) for i in range(len(gf_struct[ibl][1])) ] )

        return QN


    def fit_tails(self): 
	"""Fits the tails using the constant value for the Re Sigma calculated from F=Sigma*G.
           Works only for blocks of size one."""
	
	#if (len(self.gf_struct)==2*self.n_orb):
        if (self.blocssizeone):
            spinblocs = [v for v in self.map]
            mpi.report("Fitting tails manually")
	
            known_coeff = numpy.zeros([1,1,2],numpy.float_)
            msh = [x.imag for x in self.G[self.map[spinblocs[0]][0]].mesh ]
            fit_start = msh[self.fitting_Frequency_Start]
            fit_stop = msh[self.N_Frequencies_Accumulated]	
            
            # Fit the tail of G just to get the density
            for n,g in self.G:
                g.fitTail([[[0,0,1]]],7,fit_start,2*fit_stop) 
            densmat = self.G.density()

            for sig1 in spinblocs:
                for i in range(self.n_orb):

                    coeff = 0.0

                    for sig2 in spinblocs:
                        for j in range(self.n_orb):
                            if (sig1==sig2):
                                coeff += self.U[self.offset+i,self.offset+j] * densmat[self.map[sig1][j]][0,0].real
                            else:
                                coeff += self.Up[self.offset+i,self.offset+j] * densmat[self.map[sig2][j]][0,0].real

                    known_coeff[0,0,1] = coeff
                    self.Sigma[self.map[sig1][i]].fitTail(fixed_coef = known_coeff, order_max = 3, fit_start = fit_start, fit_stop = fit_stop)

        else:

            for n,sig in self.Sigma:

                known_coeff = numpy.zeros([sig.N1,sig.N2,1],numpy.float_)
                msh = [x.imag for x in sig.mesh]
                fit_start = msh[self.fitting_Frequency_Start]
                fit_stop  = msh[self.N_Frequencies_Accumulated]
            
                sig.fitTail(fixed_coef = known_coeff, order_max = 3, fit_start = fit_start, fit_stop = fit_stop)

		


	
def set_U_matrix(U_interact,J_hund,n_orb,l,use_matrix=True,T=None,sl_int=None,use_spinflip=False,dim_reps=None,irep=None):
    """ Set up the interaction vertex""" 

    offset = 0
    U4ind = None
    U = None
    Up = None
    if (use_matrix):
        if not (sl_int is None):
            Umat = Umatrix(l=l)
            assert len(sl_int)==(l+1),"sl_int has the wrong length"
            if (type(sl_int)==ListType):
                Rcl = numpy.array(sl_int)
            else:
                Rcl = sl_int
            Umat(T=T,Rcl=Rcl)
        else:
            if ((U_interact==None)and(J_hund==None)):
                mpi.report("Give U,J or Slater integrals!!!")
                assert 0
            Umat = Umatrix(U_interact=U_interact, J_hund=J_hund, l=l)
            Umat(T=T)
            
        Umat.reduce_matrix()
        if (Umat.N==Umat.Nmat):
            # Transformation T is of size 2l+1
            U = Umat.U
            Up = Umat.Up
        else:
            # Transformation is of size 2(2l+1)
            U = Umat.U
         # now we have the reduced matrices U and Up, we need it for tail fitting anyways

        if (use_spinflip):
            #Take the 4index Umatrix
            # check for imaginary matrix elements:
            if (abs(Umat.Ufull.imag)>0.0001).any():
                mpi.report("WARNING: complex interaction matrix!! Ignoring imaginary part for the moment!")
                mpi.report("If you want to change this, look into Wien2k/solver_multiband.py")
            U4ind = Umat.Ufull.real
    
        # this will be changed for arbitrary irep:
        # use only one subgroup of orbitals?
        if not (irep is None):
            #print irep, dim_reps
            assert not (dim_reps is None), "Dimensions of the representatives are missing!"
            assert n_orb==dim_reps[irep-1],"Dimensions of dimrep and n_orb do not fit!"
            for ii in range(irep-1):
                offset += dim_reps[ii]
    else:
        if ((U_interact==None)and(J_hund==None)):
            mpi.report("For Kanamori representation, give U and J!!")
            assert 0
        U  = numpy.zeros([n_orb,n_orb],numpy.float_)
        Up = numpy.zeros([n_orb,n_orb],numpy.float_)
        for i in range(n_orb):
            for j in range(n_orb):
	        if (i==j):
	            Up[i,i] = U_interact + 2.0*J_hund
	        else:
	       	    Up[i,j] = U_interact
		    U[i,j]  = U_interact - J_hund

    return U, Up, U4ind, offset
