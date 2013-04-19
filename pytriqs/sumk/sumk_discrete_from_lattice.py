
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

from sumk_discrete import SumkDiscrete
from pytriqs.lattice.tight_binding import TBLattice

class SumkDiscreteFromLattice (SumkDiscrete) :
    r"""
      * Computes 
      
       .. math::
         G \leftarrow \sum_k (\omega + \mu - \epsilon_k - \Sigma(k,\omega))^{-1} 
     
       for GF functions with blocks of the size of the matrix eps_k with a discrete sum.
      
      * The object contains the discretized hoppings and points in the arrays
        Hopping, BZ_Points,BZ_weights,Mu_Pattern,Overlap (IF non orthogonal)
        It can also generate a grid (ReComputeGrid) for a regular grid or a Gauss-Legendre sum
        for the whole Brillouin Zone or a patch of the BZ.
    """
    
    def __init__(self, lattice, patch = None, n_points = 8, method = "Riemann") : 
	"""   
        :param lattice: The underlying pytriqs.lattice or pytriqs.super_lattice provinding t(k)
        :param n_points:  Number of points in the BZ in EACH direction
        :param method: Riemann (default) or 'Gauss' (not checked) 
        """
        assert isinstance(lattice,TBLattice), "lattice must be a TBLattice instance"
        self.SL = lattice
        self.patch,self.method = patch,method
        # init the array
        SumkDiscrete.__init__ (self, dim = self.SL.dim, gf_struct = lattice.OrbitalNames)
        self.Recompute_Grid(n_points,  method)

     #-------------------------------------------------------------

    def __reduce__(self) : 
        return self.__class__,  (self.SL, self.patch, self.BZ_weights.shape[0],self.method) 

     #-------------------------------------------------------------

    def Recompute_Grid (self, n_points, method="Riemann", Q=None) :
        """(Re)Computes the grid on the patch given at construction : 
        
        * n_points :  Number of points in the BZ in EACH direction
        * method  : Riemann (default) or 'Gauss' (not checked) 
        * Q : anything from which a 1d-array can be computed.
              computes t(k+Q) instead of t(k) (useful for bare chi_0)
        """
        assert method in ["Riemann","Gauss"], "method %s is not recognized"%method
        self.method = method
        self.resize_arrays(n_points)
        if self.patch : 
            self.__Compute_Grid_One_patch(self.patch, n_points , method, Q)
        else : 
            self.__Compute_Grid(n_points, method, Q)

     #-------------------------------------------------------------

    def __Compute_Grid (self, n_bz,  method="Riemann", Q=None) :
	"""
        Internal
	"""

	n_bz_A,n_bz_B, n_bz_C = n_bz, (n_bz if self.dim > 1 else 1), (n_bz if self.dim > 2 else 1)
	nk = n_bz_A* n_bz_B* n_bz_C
        self.resize_arrays(nk)

	# compute the points where to evaluate the function in the BZ and with the weights
	def pts1d(N): 
	    for n in range(N) :
		yield (n - N/2 +1.0) / N

	if method=="Riemann" : 
	    BZ_weights=1.0/nk
	    k_index =0
	    for nz in pts1d(n_bz_C) : 
		for ny in pts1d(n_bz_B) :
		    for nx in pts1d(n_bz_A) : 
			self.BZ_Points[k_index,:] = (nx,ny,nz)[0:self.dim]
			k_index +=1

	elif method=="Gauss" : 
	    assert 0, "Gauss : NR gauleg not checked"
	    k_index =0
	    for wa,ptsa in NR.Gauleg(-pi,pi,n_bz_A) :
		for wb,ptsb in NR.Gauleg(-pi,pi,n_bz_B) :
		    for wc,ptsc in NR.Gauleg(-pi,pi,n_bz_C) :
			self.BZ_Points[k_index,:] = (ptsa,ptsb,ptsc)[0:self.dim] /(2*pi) 
			self.BZ_weights[k_index] = wa * wb * wc
			k_index +=1
	else :
	    raise IndexError, "Summation method unknown"

	# A shift 
	if Q :
            try : 
                Q = numpy.array(Q)
                assert len(Q.shape) ==1
            except :
                raise RuntimeError, "Q is not of correct type"
	    for k_index in range(self.N_kpts()) : 
		self.BZ_Points[k_index,:] +=Q

	# Compute the discretized hoppings from the Superlattice
	self.Hopping[:,:,:] = self.SL.hopping(self.BZ_Points.transpose()).transpose(2,0,1)

	if self.orthogonal_basis: 
            self.Mu_Pattern[:,:] =  self.SL.MuPattern[:,:]
	else :
	    assert 0 , "not checked"
	    self.Overlap[:,:,:] = self.SL.Overlap(BZ_Points.transpose())
	    mupat = self.SL.Mu_Pattern()
	    for k_index in range(self.N_kpts()) : 
		self.Mu_Pattern[:,:,k_index] = Num.dot( mupat ,self.Overlap[:,:,k_index])

    #-------------------------------------------------------------


    def __Compute_Grid_One_patch(self, patch, n_bz, method = "Riemann", Q=None) :
	"""
	Internal
	"""

        tritemp = numpy.array(patch._triangles)
        ntri = len(tritemp)/3
        nk = n_bz*n_bz*ntri
        self.resize_arrays(nk)

	# Reshape the list to group 3 points together
        triangles = tritemp.reshape((ntri,3,2))
        total_weight = 0

	# Loop over all k-points in the triangles
        k_index = 0
        for (a,b,c),w in zip(triangles,patch._weights):
          g = ((a+b+c)/3.0-a)/n_bz;
          for i in range(n_bz):
            s = i/float(n_bz)
            for j in range(n_bz-i):
              t = j/float(n_bz)
              for k in range(2):
                rv = a+s*(b-a)+t*(c-a)+(k+1)*g
                if k == 0 or j < n_bz-i-1:
	          self.BZ_Points[k_index]  = rv
	          self.BZ_weights[k_index] = w/(n_bz*n_bz)
                  total_weight += self.BZ_weights[k_index]
                  k_index = k_index+1

        # Normalize weights so that they sum up to 1
        self.BZ_weights /= total_weight

	# Compute the discretized hoppings from the Superlattice
	self.Hopping[:,:,:] = self.SL.hopping(self.BZ_Points.transpose()).transpose(2,0,1)

	if self.orthogonal_basis: 
            self.Mu_Pattern[:,:] =  self.SL.MuPattern[:,:]
	else :
	    assert 0 , "not checked"
