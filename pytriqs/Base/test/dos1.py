
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

from pytriqs.Base.Lattice.TightBinding import *

# Define the Bravais Lattice : a square lattice in 2d
BL = bravais_lattice(Units = [(1,0,0) , (0,1,0) ], Orbital_Positions= {"" :  (0,0,0)} ) 

# Prepare a nearest neighbour hopping on BL
t   = -1.00                # First neighbour Hopping
tp  =  0.0*t               # Second neighbour Hopping

# Hopping[ Displacement on the lattice] = [[t11,t12,t13....],[t21,t22,t23....],...,[....,tnn]] 
# where n=Number_Orbitals
hop= {  (1,0)  :  [[ t]],       
        (-1,0) :  [[ t]],     
        (0,1)  :  [[ t]],
        (0,-1) :  [[ t]],
        (1,1)  :  [[ tp]],
        (-1,-1):  [[ tp]],
        (1,-1) :  [[ tp]],
        (-1,1) :  [[ tp]]}

TB = tight_binding ( BL, hop)

# Compute the density of states
d = dos (TB, nkpts= 500, neps = 101, Name = 'dos')[0]

from pytriqs.Base.Archive import HDF_Archive
R = HDF_Archive('dos1.output.h5','w')
R['SquareLatt'] = d

