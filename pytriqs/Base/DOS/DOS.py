
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

import types,string,itertools
from operator import isSequenceType
import numpy

class DOS :
    r"""
* Stores a density of state of fermions  

.. math::
   :center: 

 \rho (\epsilon) \equiv \sum'_k \delta( \epsilon - \epsilon_k)

* The sum is normalized 

.. math::

  \int_{-\infty}^{\infty} d\epsilon \rho (\epsilon) = 1

* Implement :ref:`Plot Protocol <plotting>`.

    """
    def __init__(self, eps, rho ,Name = ''):
        """  
Parameters 
------------
eps : 1d array-type
    eps[i] is value of epsilon.
rho : 1d array-type
    The corresponding value of the dos. 
Name : string
     Name of the dos/orbital

        """
        self.Name = Name
        try :
            self.eps = numpy.array( eps )
            assert  len(self.eps.shape) ==1
        except :
            raise RuntimeError, "Argument eps mismatch"
        try :
            self.rho = numpy.array( rho )
            assert  len(self.rho.shape) ==1
        except :
            raise RuntimeError, "Argument rho mismatch"
        assert self.eps.shape[0] == self.rho.shape[0], "Dimensions of eps and rho do not match"
        
        
    #-------------------------------------------------------------

    def __reduce__(self) : 
        return self.__class__, (self.eps,self.rho, self.Name)

    def __reduce_to_dict__(self) :
        return {'epsilon' : self.eps, 'rho': self.rho, 'Name' : self.Name}

    @classmethod
    def __factory_from_dict__(cls,D) :
        return cls(D['epsilon'],D['rho'], D['Name'])
 
    def __repr__(self) : 
        return  """
        DOS object :
        """%self.__dict__ 

    def _plot_(self, Options) : 
        return  [ {'Type' : "XY", 'label' : self.Name, 'xlabel' :r'$\epsilon$','ylabel' : r'%s$(\epsilon)$'%self.Name, 'xdata' : self.eps,'ydata' : self.rho } ]

    def density(self,mu=0):
        """Calculates the density of free fermions for the given DOS for chemical potential mu."""

        dens = 0.0
        a = [ (e>mu) for e in self.eps ]
        try:
            ind = a.index(True)
        except:
            ind = self.eps.shape[0]

        de = self.eps[1]-self.eps[0]
        #for e,r in itertools.izip(self.eps[0:ind],self.rho[0:ind]):
        #    dens += r
        dens = (sum(self.rho[0:ind]) - self.rho[0]/2.0 - self.rho[ind-1]/2.0) * de
        #dens2 = dens + (self.rho[ind-1]/2.0 + self.rho[ind]/2.0) * de
        if (ind<self.eps.shape[0]): dens += (mu-self.eps[ind-1]) * (self.rho[ind-1] + self.rho[ind])/2.0
        return dens 

##########################################################################

def DOS_from_file(Filename, Name = '', OneOrbitalOnly = None):
    """   
    Read the DOS from a file 

    :param Filename:  a string  : name of the file
    :param Name: name of the DOS
    :param OneOrbitalOnly: can be None or an integer.
                    
    :rtype: 
       * if OneOrbitalOnly== None, returns a tuple of DOS (even if there is one dos !).
       * If OneOrbitalOnly==i, return only ONE DOS corresponding to ith orbital (starting at 1).

    Format of the file :   
        * N_orbitals +1 columns, 
        * the first column is the value of epsilon
        * the N_orbitals other columns are the values of the dos for various orbitals
    """
    f = open(Filename); s=''
    while not(s.strip()) :
        s= f.readline()
        assert s, "File is empty !"
    N_Orbitals = len (s.split()) - 1
    assert N_Orbitals >0, "File : wrong format"
    # not very safe :  fromfile routine can crashes if given non numerics
    r = numpy.fromfile(Filename,sep=' ')
    l,div  =  r.shape[0], N_Orbitals +1 
    assert l%(div)==0,"File does not contains N*%d numbers !"%(div)
    r.shape =  l//(div) , div # reshape the array
    eps = r[:,0]
    if OneOrbitalOnly : 
        assert OneOrbitalOnly>0 and OneOrbitalOnly <= N_Orbitals, " OneOrbitalOnly  "
        return DOS (r[:,0] ,r[:,OneOrbitalOnly], Name)
    else :
        return [  DOS (r[:,0] ,r[:,i +1 ], Name) for i in range (N_Orbitals)]


##########################################################################

class DOS_from_function(DOS):
    """
    * A DOS class, but constructed from a function.
    
    * The number of points can be variable and self-adjusted in the Hilbert transform to adapt precision.
    """
    def __init__(self, Function, xmin,xmax, Npts= 100,Name = ''):
        """
        :param Function: * a function :math:`\\epsilon \\rightarrow \\rho(\\epsilon)`
                         * The result type can be anything from which a 1d-array can be constructed by numpy
        :param xmin,xmax: Bound of the mesh (domain of the function). 
        :param Npts: Number of points in the mesh.
        :param Name: Name of the DOS.
        """
        assert callable(Function), "Function is not callable"
        self.Function,self.xmin,self.xmax = Function,xmin,xmax
        try :
            e = Function(0.001)
            len(numpy.array(e).shape) ==1
        except :
            raise RuntimeError , "Value of the function must be a 1d-array"
        self.__f(Npts) # compute arrays
        DOS.__init__(self,self.eps,self.rho,Name) 
        
    #-------------------------------------------------------------
    
    def __reduce__(self) : 
        return  self.__class__, (self.Function,self.xmin, self.xmax, len(self.eps), self.Name)
    
    #-------------------------------------------------------------
  
    def __f(self,N) :
        r = (self.xmax - self.xmin)/float(N-1)
        self.eps  = numpy.array( [self.xmin + r* i for i in range(N) ] )
        self.rho  = numpy.array( [self.Function(e) for e in self.eps])

#-----------------------------------------------------
#  Register the class for HDF_Archive
#-----------------------------------------------------

from pytriqs.Base.Archive.HDF_Archive_Schemes import register_class
register_class (DOS)

