
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

__all__=['TBSuperLattice']
import numpy
from pytriqs.base.Utility.myUtils import nint_strict

from TightBinding import TBLattice as Lattice

class TBSuperLattice (Lattice):
    """
    Builds a superlattice on top of a base Lattice, by picking superlattice units and optionally the cluster points.
    """
    def __init__(self,baseLattice,SuperLatticeUnits, ClusterSites = None, Remove_Internal_Hopping_in_Cluster = False ): 
        """
         * baseLattice : The underlying lattice

         * SuperLatticeUnits : The unit vectors of the superlattice in the baseLattice (integer) coordinates

	 * ClusterSites [None] : Coordinates of the cluster in baseLattice coordinates.
                   If None, an automatic computation of cluster positions is made as follows : 
                   it takes all points whose coordinates in the basis of the superlattice are in [0,1[^dimension

         * Remove_Internal_Hopping_in_Cluster : If true, the hopping terms are removed inside the cluster. 
                   Useful to add Hartree Fock terms at the boundary of a cluster, e.g.
        """
        if not isinstance(baseLattice,Lattice) : raise ValueError, "BaseLattice should be an instance of Lattice"
        self.__baseLattice = BaseLattice
        dim = baseLattice.dim

        try : 
            self.__SuperLatticeUnits = numpy.array(SuperLatticeUnits, copy=True)
            assert self.__SuperLatticeUnits.shape == (dim,dim)
        except : 
            raise ValueError,"SuperLatticeUnits is not correct. Cf Doc"

        NClusterSites = int(numpy.rint(abs(numpy.linalg.det(self.__SuperLatticeUnits ))))
        assert NClusterSites >0,"Superlattice vectors are not independant !"
        self._M = self.__SuperLatticeUnits.transpose()
        self._Mtilde = numpy.array(numpy.rint(numpy.linalg.inv(self._M)*NClusterSites), dtype = int)

        self.__Remove_Internal_Hopping_in_Cluster = Remove_Internal_Hopping_in_Cluster
        #self.norb = baseLattice.NOrbitalsInUnitCell
        self.Norb = baseLattice.NOrbitalsInUnitCell * NClusterSites
        
        # ClusterSites computation
        if ClusterSites!=None:
            self.__ClusterSites = list(ClusterSites)[:]
        else: # Computes the position of the cluster automatically 
            self.__ClusterSites = []
            #We tile the super-cell with the baseLattice points and retains
            # the points inside it and store it.
            M=numpy.array(self.__SuperLatticeUnits)
            assert M.shape==tuple(2*[dim]),"Tiling Construct : SuperLatticeUnits does not have the correct size"
            #Minv = NClusterSites*numpy.linalg.inverse(M)  #+0.5  # +0.5 is for the round
            #Mtilde = Minv.astype(numpy.Int)  # now is integer.
            Mtilde = nint_strict(NClusterSites*numpy.linalg.inv(M))
            # round to the closest integer, with assert that precision is <1.e-9
            if dim==1   :  a=(max(M[0,:]),0,0 )
            elif dim==2 :  a=(2*max(M[0,:]),2*max(M[1,:]),0 )
            elif dim==3 : a= (3*max(M[0,:]),3*max(M[1,:]), 3*max(M[2,:]))
            else : raise ValueError, "dim is not between 1 and 3 !!"
            r = lambda i :  range(-a[i] , a[i]+1)
            for nx in r(0):
                for ny in r(1):
                    for nz in r(2):
                        nv = numpy.array([nx,ny,nz][0:dim])
                        k_sl = numpy.dot(Mtilde,nv)
                        if (min(k_sl) >= 0) and (max(k_sl) < NClusterSites ) : # The point is a point of the cluster. We store it. 
                            self.__ClusterSites.append(nv.tolist())
                            
        assert len(self.__ClusterSites) == NClusterSites, """Number of cluster positions incorrect (compared to the volume of unit cell of the Superlattice)"""
        self.NClusterSites = NClusterSites

        # creating a dictionnary position_of_sites -> number e.g. (1,0) : 2 etc...
        # self._clustersites_hash =  dict ([ (tuple(pos),n) for n,pos in enumerate(self.ClusterSites)])

        # Compute the new Hopping in the supercell
        Hopping = self.Fold(baseLattice.HoppingDictionnary(), Remove_Internal_Hopping_in_Cluster)
        if 0 : 
            for k,v in Hopping.items() :
                print k
                print v.real

        # Compute the new units of the lattice in real coordinates
        Units = numpy.dot(self.__SuperLatticeUnits,baseLattice.Units)

        # Positions and names of orbitals in the supercell : just translate all orbitals for cluster site positions
        # in R^3 coordinates.
        Orbital_Positions = [POS + baseLattice.latt_to_real_x(CS) for POS in BaseLattice.OrbitalPositions for CS in self.__ClusterSites]

        #Orbital_Names = [ '%s%s'%(n,s) for n in baseLattice.OrbitalNames for s in range(NClusterSites)]
        site_index_list, orbital_index_list = range(1,NClusterSites+1),baseLattice.OrbitalNames
        if len(orbital_index_list)==1 :
            Orbital_Names= [ s for s in site_index_list ] 
        elif len(site_index_list)==1 and len(orbital_index_list)>1 :
            Orbital_Names= [ o for o in orbital_index_list]
        elif len(site_index_list)>1 and len(orbital_index_list)>1 :
            Orbital_Names= [ (pos,o) for o in orbital_index_list for pos in site_index_list]

        #print baseLattice.OrbitalNames #Orbital_Names
        Lattice.__init__(self,Units,Hopping, Orbital_Positions, Orbital_Names)
        # we pass False since the folding has arealdy been done in baseLattice

        assert self.Norb == self.NOrbitalsInUnitCell

    __HDF_reduction__ = ['__baseLattice', '__SuperLatticeUnits', '__ClusterSites','__Remove_Internal_Hopping_in_Cluster']

    def __reduce__ (self) : 
        return tuple([getattr(self,x) for x in self.__HDF_reduction__])

    def Fold(self, D1, remove_internal=False, CreateZero = None):
        """ Input : a function  r-> f(r) on the baseLattice given as a dictionnary
            Output : the function R-> F(R) folded on the superlattice.
            Only requirement is that f(r)[orbital1,orbital2] is properly defined.
            Hence f(r) can be a numpy, a GFBloc, etc...
            """
        #Res , norb = {} , self.__baseLattice.NOrbitalsInUnitCell 
        Res , norb = {} , len(D1.values()[0])
        pack = self.pack_IndexSiteOrbital
        for nsite,CS in enumerate(self.__ClusterSites) : 
            for disp,t in D1.items():
                R, alpha = self.ChangeCoordinates_L_to_SL(numpy.array(CS)+numpy.array(disp))
                if R not in Res : Res[R] = CreateZero() if CreateZero else numpy.zeros( (self.Norb,self.Norb), dtype = numpy.complex_)
                if not(remove_internal) or R!= self.baseLattice.dim*(0,): 
                    for orb1 in range(norb): 
                        for orb2 in range(norb):
                            Res[R][pack(nsite,orb1), pack(alpha,orb2)] += t[orb1,orb2]
        return Res

    def ChangeCoordinates_SL_to_L ( self, R , alpha) :
        """Given a point in the supercell R, site (number) alpha, it computes its position on the baseLattice in lattice coordinates"""
        return numpy.dot (self._M,numpy.array(R)) + self.__ClusterSites[alpha,:]
    
    def ChangeCoordinates_L_to_SL  (self, x) : 
        """Given a point on the baseLattice in lattice coordinates, returns its coordinates (R,alpha) in the Superlattice"""
        aux  = numpy.dot(self._Mtilde,numpy.array(x))
        R = aux // self.NClusterSites 
        dx = list (x - numpy.dot (self._M,R) ) # force int ?
        return tuple(R), self.__ClusterSites.index(dx) 
        
    def pack_IndexSiteOrbital(self, n_site, n_orbital) :
        """ nsite and n_orbital must start at 0"""
        return n_site + (n_orbital ) * self.NClusterSites

    def unpack_IndexSiteOrbital (self,index):
        """Inverse of pack_IndexSiteOrbital"""
	n_orbital  =   (index)//self.NClusterSites
	n_site =  index - n_orbital*self.NClusterSites
	return n_site, n_orbital

    def ClusterSites(self):
        """
           Generate the position of the cluster site in the baseLattice coordinates.
        """
        for pos in self.__ClusterSites : 
            yield pos 

    def __repr__(self) : 
        def f(A) :
            return list([ tuple(x) for x in A])
        return """SuperLattice class : \n
   base Lattice : %s
   SuperLattice Units : %s
   Remove_Internal_Hopping_in_Cluster : %s
   Cluster Sites Positions %s"""%(self.__baseLattice,f(self.__SuperLatticeUnits),self.__ClusterSites,self.__Remove_Internal_Hopping_in_Cluster)
    
 
