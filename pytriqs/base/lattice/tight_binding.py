
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

__all__ = [ 'bravais_lattice','tight_binding','dos','dos_patch','energies_on_bz_grid','energies_on_bz_path','hopping_stack']

from pytriqs_LatticeTools import bravais_lattice, tight_binding,dos_patch as dos_patch_c, dos as dos_c, energies_on_bz_grid,energies_on_bz_path,hopping_stack
from pytriqs.base.dos import DOS

def dos( TB, nkpts, neps, name) : 
    """
    :param TB: a tight_binding object
    :param nkpts: the number of k points to use in each dimension
    :param neps: number of points used in the binning of the energy
    :param name: name of the resulting dos

    :rtype: return a list of DOS, one for each band
    """
    eps, arr = dos_c(TB, nkpts,neps)
    return [ DOS (eps, arr[:,i], name) for i in range (arr.shape[1]) ]

def dos_patch( TB, triangles, nkpts, ndiv, name) :  
    """
    To be written
    """
    eps, arr = dos_c(TB, nkpts,eps)
    return DOS (eps, arr, name)




import numpy

# for backward compatibility. Not documented. 
class TBLattice : 

    def __init__ (self, Units, Hopping, Orbital_Positions= [ (0,0,0) ] , Orbital_Names = [""]) : 

        # the k are int32 which boost python does like to convert 
        def reg(k) : return tuple( int(x) for x in k) 
        self._hop = dict ( ( reg(k),numpy.array(v)) for k,v in Hopping.items())
        orb = dict ( (str(i),orb) for (i,orb) in enumerate(Orbital_Positions))
        self.bl = bravais_lattice(Units, orb)
        self.tb = tight_binding ( self.bl, self._hop) #, Orbital_Positions)
        self.dim = self.bl.dim
        self.NOrbitalsInUnitCell = self.bl.n_orbitals
        self.Units =Units
        self.OrbitalPositions = Orbital_Positions
        self.OrbitalNames = Orbital_Names
        self.MuPattern = numpy.identity(self.NOrbitalsInUnitCell)

    def latt_to_real_x(self, p) : 
        return self.bl.lattice_to_real_coordinates (numpy.array(p, numpy.float64))
        # modified since array are not converted automatically any more
        ##return self.bl.lattice_to_real_coordinates (p ) #numpy.array(p.float64))

    def HoppingDictionnary(self) : return self._hop

    def Hopping(self,k_stack) :
        return hopping_stack(self.tb,k_stack)

    #def dos(self) : d = dos (TB, nkpts= 100, neps = 100, name = 'dos2')


