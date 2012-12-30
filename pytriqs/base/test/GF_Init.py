
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
from pytriqs.base.gf_local.GFBloc_ReFreq import *
from pytriqs.base.gf_local.GFBloc_ImFreq import *
from pytriqs.base.gf_local.GF import GF
import pytriqs.base.gf_local.GF_Initializers as GF_Initializers
import numpy



h=HDF_Archive('GF_Init.output.h5','w')

g = GFBloc_ImFreq(Indices = ['eg1','eg2'], Beta = 50, NFreqMatsubara = 100, Name = "egBlock")
g['eg1','eg1'] <<= GF_Initializers.SemiCircular(HalfBandwidth = 1)
g['eg2','eg2'] <<= GF_Initializers.SemiCircular(HalfBandwidth = 2)


h['g1'] = g

g <<= GF_Initializers.Const(numpy.array([[1,2],[2,3]]))

h['g2'] = g

some_mesh = numpy.arange(-5,5,0.1)
g = GFBloc_ReFreq(Indices = ['eg1','eg2'], Beta = 50, MeshArray = some_mesh, Name = "egBlock")

g['eg1','eg1'] <<= GF_Initializers.A_Omega_Plus_B(1.0,-1.0)
g['eg2','eg2'] <<= GF_Initializers.A_Omega_Plus_B(1.0,1.0)

h['g3'] = g

del h
