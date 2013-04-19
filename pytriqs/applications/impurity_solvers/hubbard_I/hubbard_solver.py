
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

from types import *
from pytriqs.applications.dft.U_matrix import Umatrix
from pytriqs.gf.local import *
from hubbard_I import gf_hi_fullu, sigma_atomic_fullu
import pytriqs.utility.mpi as mpi
from itertools import izip
import numpy
import copy

class Solver:
    """
       Hubbard I Solver
    """

    # initialisation:
    def __init__(self, beta, l, n_msb=1025, use_spin_orbit=False):

        self.name = "Hubbard I"
        self.beta = beta
        self.l = l
        self.Nmsb = n_msb
        self.UseSpinOrbit = use_spin_orbit
        self.Converged = False
        self.Nspin = 2

        self.Nlm = 2*l+1
        if (use_spin_orbit):
            # no blocks!
            self.gf_struct = [ ('ud', range(2*self.Nlm)) ]
        else:
            # up/down blocks:
            self.gf_struct = [ ('up', range(self.Nlm)), ('down', range(self.Nlm)) ]

        # construct Greens functions:
        self.a_list = [a for a,al in self.gf_struct]
        glist = lambda : [ GfImFreq(indices = al, beta = self.beta, n_points = self.Nmsb) for a,al in self.gf_struct]
        self.G = BlockGf(name_list = self.a_list, block_list = glist(),make_copies=False)
        self.G_Old = self.G.copy()
        self.G0 = self.G.copy()
        self.Sigma = self.G.copy()
        self.Sigma_Old = self.G.copy()

        # prepare self.ealmat
        self.ealmat = numpy.zeros([self.Nlm*self.Nspin,self.Nlm*self.Nspin],numpy.complex_)

        # Define Atomic Levels Dictionary according to the GF Bloc Structure
        self.Eff_Atomic_Levels = {}
        for a,al in self.gf_struct:
            if (self.UseSpinOrbit):
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm*2,self.Nlm*2],numpy.complex_)
            else:
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm,self.Nlm],numpy.complex_)



    def solve(self, U_int, J_hund, T=None, verbosity=0, Iteration_Number=1, Test_Convergence=0.0001):
        """Calculation of the impurity Greens function using Hubbard-I"""

        if self.Converged :
            mpi.report("Solver %(name)s has already converged: SKIPPING"%self.__dict__)
            return

        self.verbosity = verbosity
        self.Nmoments = 5

        # U matrix:
        # l = (Nlm-1)/2
        # If T is specified, it is used to transform the Basis set
        Umat = Umatrix(U_interact=U_int, J_hund=J_hund, l=self.l)
        Umat(T=T)
        Umat.reduce_matrix()
        assert (Umat.N==Umat.Nmat),"Transformation that mixes spins is not implemented in hubbard_I Solver!!"

        # now we have the reduced matrices U and Up
        self.ur = Umat.Ufull
        self.umn  = Umat.Up             # reduced matrix, opposite spins
        self.ujmn = Umat.U              # reduced matrix, same spins

        M = [x for x in self.G.mesh]
        self.zmsb = numpy.array([x for x in M],numpy.complex_)

        # for the tails:
        self.tailtempl={}
        for sig,g in self.G:
            self.tailtempl[sig] = copy.deepcopy(g.tail)
            for i in range(9): self.tailtempl[sig][i] *= 0.0

        self.__save_eal('eal.dat',Iteration_Number)

        mpi.report( "Starting Fortran solver %(name)s"%self.__dict__)

        self.Sigma_Old <<= self.Sigma
        self.G_Old <<= self.G

        # call the fortran solver:
        temp = 1.0/self.beta
        gf,tail,self.atocc,self.atmag = gf_hi_fullu(e0f=self.ealmat, ur=self.ur, umn=self.umn, ujmn=self.ujmn,
                                                    zmsb=self.zmsb, nmom=self.Nmoments, ns=self.Nspin, temp=temp, verbosity = self.verbosity)

        #self.sig = sigma_atomic_fullu(gf=self.gf,e0f=self.eal,zmsb=self.zmsb,ns=self.Nspin,nlm=self.Nlm)

        if (self.verbosity==0):
            # No fortran output, so give basic results here
            mpi.report("Atomic occupancy in Hubbard I Solver  : %s"%self.atocc)
            mpi.report("Atomic magn. mom. in Hubbard I Solver : %s"%self.atmag)

        # transfer the data to the GF class:
        if (self.UseSpinOrbit):
            nlmtot = self.Nlm*2         # only one block in this case!
        else:
            nlmtot = self.Nlm

        M={}
        isp=-1
        for a,al in self.gf_struct:
            isp+=1
            M[a] = numpy.array(gf[isp*nlmtot:(isp+1)*nlmtot,isp*nlmtot:(isp+1)*nlmtot,:]).transpose(2,0,1).copy()
            for i in range(min(self.Nmoments,8)):
                self.tailtempl[a][i+1] = tail[i][isp*nlmtot:(isp+1)*nlmtot,isp*nlmtot:(isp+1)*nlmtot]

        glist = lambda : [ GfImFreq(indices = al, beta = self.beta, n_points = self.Nmsb, data =M[a], tail =self.tailtempl[a])
                           for a,al in self.gf_struct]
        self.G = BlockGf(name_list = self.a_list, block_list = glist(),make_copies=False)

        # Self energy:
        self.G0 <<= iOmega_n

        M = [ self.ealmat[isp*nlmtot:(isp+1)*nlmtot,isp*nlmtot:(isp+1)*nlmtot] for isp in range((2*self.Nlm)/nlmtot) ]
        self.G0 -= M
        self.Sigma <<= self.G0 - inverse(self.G)

        # invert G0
        self.G0.invert()

        def test_distance(G1,G2, dist) :
            def f(G1,G2) :
                print abs(G1.data - G2.data)
                dS = max(abs(G1.data - G2.data).flatten())
                aS = max(abs(G1.data).flatten())
                return dS <= aS*dist
            return reduce(lambda x,y : x and y, [f(g1,g2) for (i1,g1),(i2,g2) in izip(G1,G2)])

        mpi.report("\nChecking Sigma for convergence...\nUsing tolerance %s"%Test_Convergence)
        self.Converged = test_distance(self.Sigma,self.Sigma_Old,Test_Convergence)

        if self.Converged :
            mpi.report("Solver HAS CONVERGED")
        else :
            mpi.report("Solver has not yet converged")


    def GF_realomega(self, ommin, ommax, N_om, broadening=0.01):
        """Calculates the GF and spectral function on the real axis."""

        delta_om = (ommax-ommin)/(1.0*(N_om-1))

        omega = numpy.zeros([N_om],numpy.complex_)

        for i in range(N_om):
            omega[i] = ommin + delta_om * i + 1j * broadening

        temp = 1.0/self.beta
        gf,tail,self.atocc,self.atmag = gf_hi_fullu(e0f=self.ealmat, ur=self.ur, umn=self.umn, ujmn=self.ujmn,
                                                    zmsb=omega, nmom=self.Nmoments, ns=self.Nspin, temp=temp, verbosity = self.verbosity)

        for sig in self.a_list:
            for i in range(9): self.tailtempl[sig][i] *= 0.0

        # transfer the data to the GF class:
        if (self.UseSpinOrbit):
            nlmtot = self.Nlm*2         # only one block in this case!
        else:
            nlmtot = self.Nlm

        M={}
        isp=-1
        for a,al in self.gf_struct:
            isp+=1
            M[a] = numpy.array(gf[isp*nlmtot:(isp+1)*nlmtot,isp*nlmtot:(isp+1)*nlmtot,:]).transpose(2,0,1).copy()
            for i in range(min(self.Nmoments,8)):
                self.tailtempl[a][i+1] = tail[i][isp*nlmtot:(isp+1)*nlmtot,isp*nlmtot:(isp+1)*nlmtot]

        glist = lambda : [ GfReFreq(indices = al, window = (ommin, ommax), n_points = N_om, data = M[a], tail = self.tailtempl[a])
                           for a,al in self.gf_struct]       # Indices for the upfolded G
        self.G = BlockGf(name_list = self.a_list, block_list = glist(),make_copies=False)

        # Self energy:
        self.G0 = self.G.copy()
        self.Sigma = self.G.copy()
        self.G0 <<= Omega + 1j*broadening

        M = [ self.ealmat[isp*nlmtot:(isp+1)*nlmtot,isp*nlmtot:(isp+1)*nlmtot] for isp in range((2*self.Nlm)/nlmtot) ]
        self.G0 -= M
        self.Sigma <<= self.G0 - inverse(self.G)
        self.Sigma.note='ReFreq'          # This is important for the put_Sigma routine!!!

        #sigmamat = sigma_atomic_fullu(gf=gf,e0f=self.ealmat,zmsb=omega,nlm=self.Nlm,ns=self.Nspin)

        #return omega,gf,sigmamat


    def __save_eal(self,Filename,it):
        f=open(Filename,'a')
        f.write('\neff. atomic levels, Iteration %s\n'%it)
        for i in range(self.Nlm*self.Nspin):
            for j in range(self.Nlm*self.Nspin):
                f.write("%12.8f  "%self.ealmat[i,j])
            f.write("\n")
        f.close()


    def set_atomic_levels(self,eal):
        """ Helps to set correctly the variables for the atomic levels from a dictionary."""

        assert (type(eal)==DictType), "Give a dictionary to set_atomic_levels!"

        cnt = 0
        self.ealmat[:,:] *= 0.0

        for ind in eal:
            self.Eff_Atomic_Levels[ind] = copy.deepcopy(eal[ind])

            if self.UseSpinOrbit:
                for ii in range(self.Nlm*2):
                    for jj in range(self.Nlm*2):
                        self.ealmat[ii,jj] = self.Eff_Atomic_Levels[ind][ii,jj]
            else:
                for ii in range(self.Nlm):
                    for jj in range(self.Nlm):
                        self.ealmat[cnt*self.Nlm + ii,cnt*self.Nlm + jj] = self.Eff_Atomic_Levels[ind][ii,jj]

            cnt += 1
