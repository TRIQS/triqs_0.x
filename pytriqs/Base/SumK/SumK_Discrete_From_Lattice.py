
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

from SumK_Discrete import *
from pytriqs.Base.Lattice.SuperLattice import TBSuperLattice as SuperLattice, Lattice 

class SumK_Discrete_From_Lattice (SumK_Discrete) :
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
    
    def __init__(self, TheLattice, Patch = None, Number_Points_in_BZ = 8, Method = "Riemann") : 
	"""   
        :param TheLattice: The underlying pytriqs.Lattice or pytriqs.SuperLattice provinding t(k)
        :param Number_Points_in_BZ:  Number of points in the BZ in EACH direction
        :param Method: Riemann (default) or 'Gauss' (not checked) 
        """
        assert isinstance(TheLattice,Lattice), "TheLattice must be a Lattice"
        self.SL = TheLattice
        self.Patch,self.Method = Patch,Method
        # init the array
        SumK_Discrete.__init__ (self, dim = self.SL.dim, BS = TheLattice.OrbitalNames)
        self.Recompute_Grid(Number_Points_in_BZ,  Method)

     #-------------------------------------------------------------

    def __reduce__(self) : 
        return self.__class__,  (self.SL, self.Patch, self.BZ_weights.shape[0],self.Method) 

     #-------------------------------------------------------------

    def Recompute_Grid (self, Number_Points_in_BZ, Method = "Riemann", Q=None) :
        """(Re)Computes the grid on the Patch given at construction : 
        
        * Number_Points_in_BZ :  Number of points in the BZ in EACH direction
        * Method  : Riemann (default) or 'Gauss' (not checked) 
        * Q : anything from which a 1d-array can be computed.
              computes t(k+Q) instead of t(k) (useful for bare chi_0)
        """
        assert Method in ["Riemann","Gauss"], "Method %s is not recognized"%Method
        self.Method = Method
        self.resize_arrays(Number_Points_in_BZ)
        if self.Patch : 
            self.__Compute_Grid_One_Patch(self.Patch, Number_Points_in_BZ , Method, Q)
        else : 
            self.__Compute_Grid(Number_Points_in_BZ, Method, Q)

     #-------------------------------------------------------------

    def __Compute_Grid (self, N_BZ,  Method = "Riemann", Q=None) :
	"""
        Internal
	"""

	N_BZ_A,N_BZ_B, N_BZ_C = N_BZ, (N_BZ if self.dim > 1 else 1), (N_BZ if self.dim > 2 else 1)
	nk = N_BZ_A* N_BZ_B* N_BZ_C
        self.resize_arrays(nk)

	# compute the points where to evaluate the function in the BZ and with the weights
	def pts1d(N): 
	    for n in range(N) :
		yield (n - N/2 +1.0) / N

	if Method=="Riemann" : 
	    BZ_weights=1.0/nk
	    k_index =0
	    for nz in pts1d(N_BZ_C) : 
		for ny in pts1d(N_BZ_B) :
		    for nx in pts1d(N_BZ_A) : 
			self.BZ_Points[k_index,:] = (nx,ny,nz)[0:self.dim]
			k_index +=1

	elif Method=="Gauss" : 
	    assert 0, "Gauss : NR gauleg not checked"
	    k_index =0
	    for wa,ptsa in NR.Gauleg(-pi,pi,N_BZ_A) :
		for wb,ptsb in NR.Gauleg(-pi,pi,N_BZ_B) :
		    for wc,ptsc in NR.Gauleg(-pi,pi,N_BZ_C) :
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
	self.Hopping[:,:,:] = self.SL.Hopping(self.BZ_Points.transpose()).transpose(2,0,1)

	if self.Orthogonal_Basis: 
            self.Mu_Pattern[:,:] =  self.SL.MuPattern[:,:]
	else :
	    assert 0 , "not checked"
	    self.Overlap[:,:,:] = self.SL.Overlap(BZ_Points.transpose())
	    mupat = self.SL.Mu_Pattern()
	    for k_index in range(self.N_kpts()) : 
		self.Mu_Pattern[:,:,k_index] = Num.dot( mupat ,self.Overlap[:,:,k_index])

    #-------------------------------------------------------------


    def __Compute_Grid_One_Patch(self, Patch, N_BZ, Method = "Riemann", Q=None) :
	"""
	Internal
	"""

        tritemp = numpy.array(Patch._triangles)
        ntri = len(tritemp)/3
        nk = N_BZ*N_BZ*ntri
        self.resize_arrays(nk)

	# Reshape the list to group 3 points together
        triangles = tritemp.reshape((ntri,3,2))
        total_weight = 0

	# Loop over all k-points in the triangles
        k_index = 0
        for (a,b,c),w in zip(triangles,Patch._weights):
          g = ((a+b+c)/3.0-a)/N_BZ;
          for i in range(N_BZ):
            s = i/float(N_BZ)
            for j in range(N_BZ-i):
              t = j/float(N_BZ)
              for k in range(2):
                rv = a+s*(b-a)+t*(c-a)+(k+1)*g
                if k == 0 or j < N_BZ-i-1:
	          self.BZ_Points[k_index]  = rv
	          self.BZ_weights[k_index] = w/(N_BZ*N_BZ)
                  total_weight += self.BZ_weights[k_index]
                  k_index = k_index+1

        # Normalize weights so that they sum up to 1
        self.BZ_weights /= total_weight

	# Compute the discretized hoppings from the Superlattice
	self.Hopping[:,:,:] = self.SL.Hopping(self.BZ_Points.transpose()).transpose(2,0,1)

	if self.Orthogonal_Basis: 
            self.Mu_Pattern[:,:] =  self.SL.MuPattern[:,:]
	else :
	    assert 0 , "not checked"
