
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

from pytriqs.Base.GF_Local import GF_Initializers
import types,string,inspect,itertools
from operator import isSequenceType
from pytriqs.Base.DOS.DOS import *
import pytriqs.Base.Utility.MPI as MPI

class Hilbert_Transform : 
    r"""
    Computes the Hilbert Transform from a DOS object 

    .. math::
    
       \int_{-\infty}^\infty d \epsilon \rho(\epsilon) \Bigl(  (\omega + \mu +
       I\eta)\mathbf{1} - \hat\varepsilon(\epsilon) - \text{Field} - \Sigma(\epsilon)
       \Bigr)^{-1}

    """
    def __init__(self, Rho):
        """
        :param Rho: a DOS object.
        """
        self.dos  = Rho
        assert isinstance(Rho, DOS),  "See Doc. Rho must be a DOS"
        self.__normalize()

    #-------------------------------------------------------------

    def __reduce__(self) : 
        return self.__class__, (self.Rho)

    #-------------------------------------------------------------
   
    def __normalize(self):
        # normalisation. dos is not the value of the function, is the weight of the integrals
        R = numpy.array(self.dos.rho,copy=True)
        self.rho_for_sum = R
        eps = self.dos.eps
        R[0]  *= (eps[1] - eps[0])
        R[-1] *= (eps[-1] - eps[-2])
        for i in xrange(1, eps.shape[0] - 1) : 
            R[i] *=  (eps[i+1] - eps[i])/2+(eps[i] - eps[i-1])/2
        R /= numpy.sum(R) 
     
    #-------------------------------------------------------------

    def __call__ (self, Sigma, mu=0, eta = 0, Field = None, Epsilon_Hat=None, Res = None,
                  Number_Points_in_integral=None, Test_Convergence = None):
        r""" 
        Compute the Hilbert Transform 
               
        Parameters 
        -----------
        
        mu  : float
        eta : float 
        Sigma : a GFBloc or a function epsilon-> GFBloc
        Field : anything that can added to the GFBloc Sigma, e.g. : 
                 * an Array_with_GFBloc_Indices (same size as Sigma) 
                 * a GBloc 
        Epsilon_Hat : a function that takes a 1d array eps[i] and returns 3d-array   eps[i,:,:]
                            where the :,: has the matrix structure of Sigma. Default : eps[i] * Identity_Matrix
                            Used only when DOS is a DOS_from_Function : 
        Number_Points_in_integral : How many points to use. If None, use the Npts of construction
        Test_Convergence : If defined, it will refine the grid until CV is reached
                          starting from Number_Points_in_integral and multiplying by 2
        
        Returns
        --------
        
        Returns the result. If provided, use Res to compute the result locally.
        """

        # we suppose here that self.eps, self.rho_for_sum such that
        # H(z) = \sum_i self.rho_for_sum[i] * (z- self.eps[i])^-1

        # Check Sigma and Res
        assert Sigma.N1==Sigma.N2, "Sigma must be square"
        if Res :
            assert Res.N1 == Sigma.N1 and Res.N2 == Sigma.N2, "Size of Res and Sigma mismatch"
        else :
            Res = Sigma.copy()

        if not( isinstance (self.dos, DOS_from_function)):
            assert Number_Points_in_integral==None and Test_Convergence == None, " Those parameters can only be used with an Dos_from_function"
        if Field !=None : 
            try : 
                Res += Field
            except : 
                assert 0,"Field can not be added to the Green function blocks !. Cf Doc"

        def HT(Res) : 
            # First compute the eps_hat array
            eps_hat = Epsilon_Hat(self.dos.eps) if Epsilon_Hat else numpy.array( [ x* numpy.identity (Sigma.N1) for x in self.dos.eps] )
            assert eps_hat.shape[0] == self.dos.eps.shape[0],"Epsilon_Hat function behaves incorrectly"
            assert eps_hat.shape[1] == eps_hat.shape[2],"Epsilon_Hat function behaves incorrectly (result not a square matrix)"
            assert Sigma.N1 == eps_hat.shape[1], "Size of Sigma and of epsilon_hat mismatch"

            Res.zero()
            Sigma_fnt = callable(Sigma)
            if Sigma_fnt : assert len(inspect.getargspec(Sigma)[0]) ==1, "Sigma function is not of the correct type. See Documentation"

            # Perform the sum over eps[i]
            tmp,tmp2 = Res.copy(),Res.copy()
            tmp <<= GF_Initializers.A_Omega_Plus_B(1,mu + eta * 1j)
            if not(Sigma_fnt) :
                tmp -= Sigma
            if Field != None : tmp -= Field
            
            # I slice all the arrays on the node. Cf reduce operation below. 
            for d,e_h,e in  itertools.izip (*[MPI.slice_array(A) for A in [self.rho_for_sum,eps_hat,self.dos.eps]]):
                tmp2.copyFrom(tmp)
                tmp2 -= e_h
                if Sigma_fnt : tmp2 -= Sigma(e)
                tmp2.invert()
                tmp2 *= d
                Res += tmp2
            # sum the Res GF of all nodes and returns the results on all nodes...
            # Cf Boost.mpi.python, collective communicator for documentation.
            # The point is that Res is pickable, hence can be transmitted between nodes without further code...
            Res <<= MPI.all_reduce(MPI.world,Res,lambda x,y : x+y)
            MPI.barrier()
        # END of HT

        def test_distance(G1,G2, dist) :
            def f(G1,G2) : 
                dS = max(abs(G1._data.array - G2._data.array).flatten())  
                aS = max(abs(G1._data.array).flatten())
                return dS <= aS*dist
            #return reduce(lambda x,y : x and y, [f(g1,g2) for (i1,g1),(i2,g2) in izip(G1,G2)])
            return f(G1,G2) # for block function, the previous one is for GF functions

        if isinstance (self.dos, DOS_from_function): 
            
            if not(Number_Points_in_integral) : # if not defined, use the defaults given at construction of the dos
                Number_Points_in_integral=  len(self.dos.eps)
            else:
                self.dos._DOS__f(Number_Points_in_integral)
                self.__normalize()

            HT(Res)

            nloop, test = 1,0
            while Test_Convergence and nloop < 10 and (nloop == 1 or test > Test_Convergence):
                if nloop>1 :
                    self.dos._DOS__f(Number_Points_in_integral)
                    self.__normalize()

                Res_old = Res.copy()
                Res= DOS.Hilbert_Transform(self,Sigma,mu,eta,Epsilon_Hat, Res)
                test = test_distance(Res,Res_old, Test_Convergence)
                Number_Points_in_integral *=2
                
        else :  # Ordinary DOS
            HT(Res)

        return Res
