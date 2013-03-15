
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

from pytriqs.archive import *
from pytriqs.applications.impurity_solvers.hubbard_I import Solver
import numpy

S = Solver(beta = 200, l=2)

eal={}
eal['up'] = -1*numpy.identity(5)
eal['down'] = -1*numpy.identity(5)
S.set_atomic_levels(eal=eal)

S.solve(U_int = 6.0, J_hund=0.6)

ar = HDFArchive('hubbard.output.h5')
ar['G'] = S.G
ar['Sigma'] = S.Sigma
del ar

