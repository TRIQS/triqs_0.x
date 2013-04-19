
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

from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *

#
#  Example of DMFT single site solution with CTQMC
#

# set up a few parameters
Half_Bandwidth= 1.0
U = 2.5
Chemical_Potential = U/2.0
beta = 100

# Construct a CTQMC solver
from pytriqs.applications.impurity_solvers.ctqmc_hyb import Solver
S = Solver(beta = beta, gf_struct = [ ('up',[1]), ('down',[1]) ])

# init the Green function
S.G <<= SemiCircular(Half_Bandwidth)

# Impose Paramagnetism
g = 0.5*(S.G['up']+S.G['down'])
for name, bloc in S.G : bloc <<= g

# Compute G0
for sig,g0 in S.G0 :
  g0 <<= inverse( iOmega_n + Chemical_Potential - (Half_Bandwidth/2.0)**2  * S.G[sig] )

# Solve
from pytriqs.applications.impurity_solvers.operators import *
S.solve(H_local = U * N('up',1) * N('down',1),
        quantum_numbers = { 'Nup' : N('up',1), 'Ndown' : N('down',1) },
        n_cycles  = 5000,
        length_cycle = 500,
        n_warmup_cycles = 5000,
        n_time_slices_delta = 10000,
        n_time_slices_gtau = 1000,
        random_name = "",
        n_legendre = 30,
        use_segment_picture = True)

# Calculation is done. Now save a few things
# Save into the shelve
Results = HDFArchive("single_site_bethe.output.h5",'w')
Results["G"] = S.G
Results["Gl"] = S.G_Legendre

