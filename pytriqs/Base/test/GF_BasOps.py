
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

from pytriqs.Base.Archive import *
from pytriqs.Base.GF_Local import *
#import pytriqs.Base.GF_Local.GF_Initializers as GF_Initializers
import numpy, copy
from pytriqs.Base.Utility.myUtils import Conjugate

h=HDF_Archive('GF_BasOps.output.h5','w')

ga = GFBloc_ImFreq(Indices = [1,2], Beta= 50, NFreqMatsubara = 100, Name = "a1Block")
gb = GFBloc_ImFreq(Indices = [1,2], Beta= 50, NFreqMatsubara = 100, Name = "b1Block")

G = GF(NameList = ('a','b'), BlockList = (ga,gb), Copy = False)

#G <<= GF_Initializers.A_Omega_Plus_B(1.0,2.0)
G <<= iOmega_n + 2.0
h['G1'] = G

# inverse:
G <<= inverse(G)

#  Density:
dens = G.total_density()
h['dens'] = dens

# FT:
f = lambda g,L : GFBloc_ImTime(Indices= g.Indices, Beta = g.Beta, NTimeSlices=L )
gt = GF(Name_Block_Generator = [ (n,f(g,200) ) for n,g in G], Copy=False, Name='gt')
for (i,gtt) in gt : gtt.setFromInverseFourierOf(G[i])

h['gt'] = gt

# Matrix operations:
ga2 = GFBloc_ImFreq(Indices = [1,2,3], Beta= 50, NFreqMatsubara = 100, Name = "a1Block")
mat = numpy.array([[1.0,0.0,1.0],[-1.0,1.0,0.0]])

ga2.from_L_G_R(mat.transpose(),ga,mat)

h['ga2'] = ga2


# conjugate:
Gc = Conjugate(G)
h['Gc'] = Gc

# Some tail stuff:
tailtempl={}
for sig,g in G: 
    tailtempl[sig] = copy.deepcopy(g._tail)
h['tt'] = tailtempl


# tranpose
g = G['a']
gt = G['a'].transpose()
gt._data.array[0,1,3] = 100
assert g._data.array[1,0,3] == 100

del h
