
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
from pytriqs.base.gf_local.gf_refreq import *
from pytriqs.base.gf_local.gf_imfreq import *
from pytriqs.base.gf_local.block_gf import BlockGf
import pytriqs.base.gf_local.gf_init as gf_init
import numpy



h=HDFArchive('gf_init.output.h5','w')

g = GfImFreq(indices = ['eg1','eg2'], beta = 50, n_matsubara = 100, name = "egBlock")
g['eg1','eg1'] <<= gf_init.SemiCircular(half_bandwidth = 1)
g['eg2','eg2'] <<= gf_init.SemiCircular(half_bandwidth = 2)


h['g1'] = g

g <<= gf_init.Const(numpy.array([[1,2],[2,3]]))

h['g2'] = g

some_mesh = numpy.arange(-5,5,0.1)
g = GfReFreq(indices = ['eg1','eg2'], beta = 50, mesh_array = some_mesh, name = "egBlock")

g['eg1','eg1'] <<= gf_init.A_Omega_Plus_B(1.0,-1.0)
g['eg2','eg2'] <<= gf_init.A_Omega_Plus_B(1.0,1.0)

h['g3'] = g

del h
