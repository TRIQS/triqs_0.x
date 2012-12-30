
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

from pytriqs.base.Archive import HDF_Archive
from pytriqs.base.GF_Local import *
from pytriqs.base.DMFT.Loop_Generic import *
from pytriqs.base.GF_Local.Descriptors import iOmega_n,SemiCircular

#
#  Example of DMFT single site solution with CTQMC
#

# Opens the results hdf5 file
Results = HDF_Archive("SingleSiteBethe.output.h5",'w')

# set up a few parameters
Half_Bandwidth= 1.0    
U = 2.5
Chemical_Potential = U/2.0
Beta = 100

# Construct a CTQMC solver
from pytriqs.solvers.operators import *        # imports the class manipulating C, C_dagger and N = C_dagger C
from pytriqs.solvers.ctqmc_hyb import Solver   # imports the solver class
S = Solver(Beta = Beta,                                                      # inverse temperature
           GFstruct = [ ('up',[1]), ('down',[1]) ],                          # Structure of the Green function 
           H_Local = U * N('up',1) * N('down',1),                            # Local Hamiltonian 
           Quantum_Numbers = { 'Nup' : N('up',1), 'Ndown' : N('down',1) },   # Quantum Numbers (operators commuting with H_Local)
           N_Cycles  = 5000,#200000,                                         # Number of QMC cycles
           Length_Cycle = 500,
           N_Warmup_Cycles = 5000,                
           N_Time_Slices_Delta= 10000,
           N_Time_Slices_Gtau = 1000,
           Random_Generator_Name = "", 
           N_Legendre_Coeffs = 30,
           #Random_Generator_Name = "mt19937",
           Use_Segment_Picture = True)

#S.Proba_Move = 0.0
#S.Global_Moves =  [ (0.01,  lambda (a,alpha,dag) : (a, alpha, dag)),
#                   (0.01,  lambda (a,alpha,dag) : ({'up' : 'down', 'down' : 'up'}[a], alpha, dag)) ]

# init the Green function
#S.G <<= GF_Initializers.SemiCircular(Half_Bandwidth)
S.G <<= SemiCircular(Half_Bandwidth)

# Derive a DMFT loop
class myloop (DMFT_Loop_Generic) : 
      def Self_Consistency(self) :
            # Impose Paramagnetism
            g = 0.5*(S.G['up']+S.G['down'])
            for name, bloc in S.G : bloc <<= g 
            
            # Compute G0
            for sig,g0 in S.G0 :
                  #g0 <<= GF_Initializers.A_Omega_Plus_B(1,Chemical_Potential)
                  #g0 <<=  iOmega_n + Chemical_Potential 
                  #g0 -= (Half_Bandwidth/2.0)**2  * S.G[sig]
                  #g0.invert()
                  g0 <<= inverse( iOmega_n + Chemical_Potential - (Half_Bandwidth/2.0)**2  * S.G[sig] )


# instanciate and run
myloop(Solver_List = S).run(N_Loops = 1)

# Calculation is done. Now save a few things
# Save into the shelve
Results["G"] = S.G
Results["Gtau"] = S.G_tau
#Results["Sigma"] = S.Sigma
# Save into some text files...
#S.G.save("GC") 
#S.Sigma.save("Sigma")

