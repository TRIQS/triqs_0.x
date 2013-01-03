from pytriqs.base.gf_local import *
from pytriqs.solvers.operators import *
from pytriqs.solvers.ctqmc_hyb import Solver

D, V, U = 1.0, 0.2, 4.0
e_f, Beta = -U/2.0, 50

# The impurity solver
S = Solver(Beta = Beta,                             # inverse temperature
           GFstruct = [ ('up',[1]), ('down',[1]) ], # Structure of the Green's function 
           H_Local = U * N('up',1) * N('down',1),   # Local Hamiltonian 
           Quantum_Numbers = {                      # Quantum Numbers 
               'Nup' : N('up',1),                   # (operators commuting with H_Local) 
               'Ndown' : N('down',1) },          
           N_Cycles = 100000,                       # Number of QMC cycles
           Length_Cycle = 200,                      # Length of one cycle 
           N_Warmup_Iteration = 10000,              # Warmup cycles
           Use_Segment_Picture = True,              # Use the segment picture
           Global_Moves = [                         # Global move in the QMC
               (0.05, lambda (a,alpha,dag) : ( {'up':'down','down':'up'}[a],alpha,dag ) ) ], 
           )

from pytriqs.base.archive import HDFArchive
import pytriqs.base.utility.mpi as mpi

for random_name in ['mt11213b','lagged_fibonacci607']:

  # Solve using random_name as a generator
  S.Random_Generator_Name = random_name
  for spin, g0 in S.G0 :
    g0 <<= inverse( iOmega_n - e_f - V**2 * Wilson(D) ) 
  S.Solve()

  # Save the results in an hdf5 file (only on the master node)
  if mpi.is_master_node():
    Results = HDFArchive("random.h5")
    Results["G_%s"%(random_name)] = S.G
    del Results

