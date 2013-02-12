
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

from pytriqs.base.archive import *
from pytriqs.base.gf.local import gf_init, BlockGf, GfImFreq, inverse
from pytriqs.base.lattice.super_lattice import TBSuperLattice
from pytriqs.base.lattice.tight_binding import TBLattice
from pytriqs.base.sumk import *
from pytriqs.base.dmft import DMFTLoopGeneric
import pytriqs.base.utility.mpi as mpi

#
# A plaquette calculation with CDMFT, and a basic Riemann sum over BZ
#

# A few variables
Beta = 10
Density_Required = None #0.85 #None
Chemical_potential  =  2.0
Field_AF, Field_F = 0.0, 0.0

#
#  Solver
#
from pytriqs.solvers.ctqmc_hyb.models import Solver_2x2_Para_Hubbard
S = Solver_2x2_Para_Hubbard(Beta = Beta, U_interact = 4.0)
# I can change some parameter in the S. namespace, overruling the defaults
# I could even change them from iteration to iteration in the loop below.
S.N_Cycles = 3000
S.N_Warmup_Cycles = 0
S.Length_Cycle = 10
S.N_Legendre_Coeffs = 30

# Lattice and Super-Lattice
t   = -1.00             # First neighbour Hopping
tp  =  0.0*t            # Second neighbour Hopping
hop= {  (1, 0)  :  [[ t]], # Hopping[ Displacement on the lattice] = [[t11, t12, t13....], [t21, t22, t23....], ..., [...., tnn]]
        (-1, 0) :  [[ t]], # where n=Number_Orbitals
        (0, 1)  :  [[ t]],
        (0, -1) :  [[ t]],
        (1, 1)  :  [[ tp]],
        (-1, -1):  [[ tp]],
        (1, -1) :  [[ tp]],
        (-1, 1) :  [[ tp]]}

L = TBLattice ( units = [(1, 0, 0) , (0, 1, 0) ], hopping = hop)
SL = TBSuperLattice(tb_lattice =L, super_lattice_units = [ (2, 0), (0, 2) ])

# SumK function that will perform the sum over the BZ
SK = SumkDiscreteFromLattice (lattice = SL, n_points = 8, method = "Riemann")

# Defines G and Sigma with a block structure compatible with the SumK function
G= BlockGf( name_block_generator = [ (s, GfImFreq(indices = SK.GFBlocIndices, mesh = S.G.mesh)) for s in ['up', 'down'] ], make_copies = False)
Sigma = G.copy()

# Init Sigma
for n, B in S.Sigma : B <<= gf_init.Const(2.0)

# Now I write my DMFT loop...
class myloop (DMFTLoopGeneric) :
   def Self_Consistency(self) :
      S.Transform_SymmetryBasis_toRealSpace (IN= S.Sigma, OUT = Sigma) # Embedding

      # Computes sum over BZ and returns density
      F = lambda mu : SK(mu = mu, Sigma = Sigma, field = None , result = G).total_density()/4

      if Density_Required :
         self.Chemical_potential = dichotomy.dichotomy(function = F,
                                                       x_init = self.Chemical_potential, y_value =Density_Required,
                                                       precision_on_y = 0.01, delta_x=0.5,  max_loops = 100,
                                                       x_name="Chemical_Potential", y_name= "Total Density",
                                                       verbosity = 3)[0]
      else:
         mpi.report("Total density  = %.3f"%F(self.Chemical_potential))

      S.Transform_RealSpace_to_SymmetryBasis (IN = G, OUT = S.G)       # Extraction
      S.G0 = inverse(S.Sigma + inverse(S.G))                           # Finally get S.G0

# Construct an instance and run.
myloop(solver_list = S, Chemical_potential  = 2.0).run(n_loops = 1)

# Opens the results shelve
if mpi.is_master_node():
  Results = HDFArchive("cdmft_4_sites.output.h5", 'w')
  Results["G"] = S.G
