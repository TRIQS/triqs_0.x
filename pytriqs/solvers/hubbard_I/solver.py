
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
from pytriqs.solvers.hubbard_I.solver_base import SolverBaseHub
from pytriqs.base.gf_local.gf_imfreq import *
from pytriqs.base.gf_local.block_gf import BlockGf
from pytriqs.dft.U_matrix import Umatrix
import copy,numpy

class Solver(SolverBaseHub):
    """
       Hubbard I Solver
    """
   
   
    # initialisation:
    def __init__(self,beta,U_int,J_hund,l,n_msb=1025,T=None, use_spin_orbit=False, verbosity=0):
        #If T is specified, it is used to transform the Basis set

        Nlm=2*l+1
        if (use_spin_orbit):
            # no blocks!
            GFstruct = [ ('ud', range(2*Nlm)) ]
        else:
            # up/down blocks:
            GFstruct = [ ('up', range(Nlm)), ('down', range(Nlm)) ]
        
        # U matrix:
        #l = (Nlm-1)/2
        Umat = Umatrix(U_interact=U_int, J_hund=J_hund, l=l)  
        Umat(T=T)
        Umat.reduce_matrix()
        assert (Umat.N==Umat.Nmat),"Transformation that mixes spins is not implemented in hubbard_I Solver!!"
        # now we have the reduced matrices U and Up

        SolverBaseHub.__init__(self, Beta=beta, GFstruct=GFstruct, Nlm=Nlm, Nmsb = n_msb, UseSpinOrbit = use_spin_orbit, Verbosity=verbosity)
        
        self.ur = Umat.Ufull
        self.umn  = Umat.Up             # reduced matrix, opposite spins
        self.ujmn = Umat.U              # reduced matrix, same spins


        # Define Atomic Levels Dictionary according to the GF Bloc Structure
        self.Eff_Atomic_Levels = {}
        for a,al in GFstruct:
            if (self.UseSpinOrbit):
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm*2,self.Nlm*2],numpy.complex_)
            else:
                self.Eff_Atomic_Levels[a] = numpy.zeros([self.Nlm,self.Nlm],numpy.complex_)
            

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

            
