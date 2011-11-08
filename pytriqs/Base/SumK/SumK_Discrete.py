
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

__all__= ['SumK_Discrete']
from pytriqs.Base.GF_Local.GF import GF
from pytriqs.Base.GF_Local import GF_Initializers
import pytriqs.Base.Utility.MPI as MPI
from itertools import *
import inspect
import copy,numpy

class SumK_Discrete :
    """
      INTERNAL USE
      The function to compute \[ G \leftarrow \sum_k (\omega + \mu - eps_k - Sigma(k,\omega))^{-1} \]
      for GF functions with blocks of the size of the matrix eps_k with a discrete sum.
      The class contains the discretized hoppings and points in the arrays
      Hopping, BZ_Points,BZ_weights,Mu_Pattern,Overlap (IF non orthogonal)
      It can also generate a grid (ComputeGrid) for a regular grid or a Gauss-Legendre sum.
    """
    def __init__ (self,dim, BS, Orthogonal_Basis = True ):
	"""  
	Just constructs the arrays, but without initializing them
	- dim is the dimension
	- BS : Indices of the Green function
	- Orthogonal_Basis : True by default
	"""
        self.__GFBLOC_Structure = copy.deepcopy(BS)
        self.Orthogonal_Basis,self.dim = Orthogonal_Basis,dim

   #-------------------------------------------------------------

    def resize_arrays (self, nk):
	"""  
	Just constructs the arrays, but without initializing them
	- nk : total number of k points
	"""
        # constructs the arrays.
        no = len(self.__GFBLOC_Structure)
        self.Hopping    = numpy.zeros([nk,no,no],numpy.complex_)   # t(k_index,a,b)
        self.BZ_Points  = numpy.zeros([nk,self.dim],numpy.float_)      # k(k_index,:)
        self.BZ_weights = numpy.ones([nk],numpy.float_)/ float(nk) # w(k_kindex) ,  default normalisation 
        self.Mu_Pattern  =  numpy.identity(no,numpy.complex_) if self.Orthogonal_Basis else numpy.zeros([no,no,nk],numpy.complex_)
        self.Overlap = numpy.array(self.Mu_Pattern, copy=True)

   #-------------------------------------------------------------
    
    def __get_GFBloc_Structure(self) :
        """Returns the ONLY block indices accepted for the G and Sigma argument of the 
        SumK function"""
        return self.__GFBLOC_Structure

    GFBlocIndices = property(__get_GFBloc_Structure)

    #-------------------------------------------------------------

    def __call__ (self, Sigma, mu=0, eta = 0, Field = None, Res = None, SelectedBlocks = () ):
	""" 
	- Computes :
	   Res <- \[ \sum_k (\omega + \mu - Field - t(k) - Sigma(k,\omega)) \]
           if Res is None, it returns a new GF with the results.
           otherwise, Res must be a GF, in which the calculation is done, and which is then returned.
           (this allows chain calculation : SK(mu = mu,Sigma = Sigma, Res = G).total_density()
           which computes the sumK into G,  and returns the density of G.
  
        - Sigma can be a X, or a function k-> X or a function k,eps ->X where  : 
	    - k is expected to be a 1d-numpy array of size self.dim of float, 
	      containing the k vector in the basis of the RBZ  (i.e.  -0.5< k_i <0.5)
            - eps is t(k)
	    - X is anything such that X[BlockName] can be added/subtracted to a GFBloc for BlockName in SelectedBlocks.
	      e.g. X can be a GF (with at least the SelectedBlocks), or a dictionnary BlockName -> array
	      if the array has the same dimension as the GF blocks (for example to add a static Sigma).

        - Field : Any k independant  Array_with_GF_Indices to be added to the GF 

        - SelectedBlocks : The calculation is done with the SAME t(k) for all blocks. If this list is not None
	  only the blocks in this list are calculated.
	  e.g. G and Sigma have block indices 'up' and 'down'. 
	       if SelectedBlocks ==None : 'up' and 'down' are calculated
	       if SelectedBlocks == ['up'] : only 'up' is calculated. 'down' is 0.

         """
        if Field : assert isinstance(Field,Array_with_GF_Indices) , " Field must be a  Array_with_GF_Indices object. Cf Example"
        S = Sigma.View_SelectedBlocks(SelectedBlocks) if SelectedBlocks else Sigma
        Gres = Res if Res else Sigma.copy() 
        G = Gres.View_SelectedBlocks(SelectedBlocks) if SelectedBlocks else Gres

        # check input
        assert self.Orthogonal_Basis, "Local_G : must be orthogonal. non ortho cases not checked."
        assert isinstance(G,GF), "G must be a GF"
        assert list(set([ g.N1 for i,g in G])) == [self.Hopping.shape[1]],"G size and hopping size mismatch"
        assert self.BZ_weights.shape[0] == self.N_kpts(), "Internal Error"
        Sigma_Nargs = len(inspect.getargspec(Sigma)[0]) if callable (Sigma) else 0
        assert Sigma_Nargs <=2 , "Sigma function is not of the correct type. See Documentation"

        #init
        G.zero()
        #tmp,tmp2 = GF(G),GF(G)
        tmp,tmp2 = G.copy(),G.copy()
        mupat = mu * self.Mu_Pattern 
        tmp <<= GF_Initializers.A_Omega_Plus_B(A=1,B=0)
        #tmp.Set_Omega()
        ##tmp += tmp.Nblocks() * [ mupat ]
        if Field : tmp -= Field 
        if Sigma_Nargs==0: tmp -= Sigma  # substract Sigma once for all

        # Loop on k points...
        for w, k, eps_k in izip(*[MPI.slice_array(A) for A in [self.BZ_weights, self.BZ_Points, self.Hopping]]):
            tmp2 <<= tmp
            #tmp2.copy_from(tmp)
            tmp2 -= tmp2.NBlocks * [eps_k -mupat ]
            #tmp2.save("tmp2_w")
            #Sigma.save("S_w")

            if Sigma_Nargs == 1: tmp2 -= Sigma (k)
            elif Sigma_Nargs ==2: tmp2 -= Sigma (k,eps_k)
            tmp2.invert()
            tmp2 *= w
            G += tmp2
            #G.save("GG1")
            #print mu,mupat,eps_k
            #assert 0
            #print G['up'][1,1]._data
        G <<= MPI.all_reduce(MPI.world,G,lambda x,y : x+y)
        MPI.barrier()

        return Res

    #-------------------------------------------------------------

    def N_kpts(self) : 
	""" Returns the number of k points"""
	return self.BZ_Points.shape[0]
