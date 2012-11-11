################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011-2012 by M. Ferrero, O. Parcollet
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

def Delta(G0) :
    """Computes the hybridization from a Green function """
    if type(G0) not in [GFBloc_ImFreq, GFBloc_ReFreq] : 
        raise RuntimeError, "Hybridization can only be computed in frequency"""
    G0 = self if self._tail.OrderMin <=-1 else inverse(self)
    tmp = G0.copy()
    tmp <<= GF_Initializers.A_Omega_Plus_B(G0._tail[-1], G0._tail[0])
    tmp -= G0
    return tmp


