
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

r"""
GFInitializer functions are used to initialize the GFBloc_xxx objects
  Usage : 
    any_GFInitializer(G, parameters)

  Action :
    Replaces the _data and/or the _tail in G (depending of the initialiser)

  Possible Initializers : 
    * Basic : simple replacement of data and tail
    * Function : from a python function : omega -> G[omega] plus the tail
    * FreeBath: a free bath with a given number of sites
    * SemiCircular :...

  Motivations : 
    Put the init out of the Base class. Hence, one can add a new initialiser
    without touching the GFBloc_xxx classes

""" 

import numpy
from math import *
from pytriqs_GF import GF_Statistic,GF_Type,TailGf,MeshGf
from pytriqs.base.utility.my_utils import sign
from pytriqs.base.gf_local.array_view import ArrayViewWithIndexConverter

class Base:
    def __init__(self,**kargs) :
        self.__dict__.update(kargs)

#########################################################################

class Function (Base):
    def __init__ (self, F, Tail = None) : 
        Base.__init__(self,F=F, Tail=Tail)
        
    def __call__(self,G) :
        if not(callable(self.F)) : raise RuntimeError, "GFInitializer.Function : f must be callable"
        res = G._data.array[:,:,:]
        try : 
            for n,om in enumerate(G.mesh) : res[:,:,n] = self.F(om)
        except :
            raise RuntimeError, "The given function has a problem..."
        if self.Tail : G._tail.copy_from(self.Tail)
        return G

#########################################################################

class Const(Base):
    def __init__ (self, C) :
        if isinstance(C,ArrayViewWithIndexConverter) : C=C.array # trasnform to numpy.array
        Base.__init__(self, C=C)
         
    def __call__(self,G) :
        C = self.C
        if G.mesh.TypeGF not in [GF_Type.Imaginary_Frequency, GF_Type.Real_Frequency] : 
            raise TypeError, "This initializer is only correct in frequency"

        if not isinstance(C,numpy.ndarray) : 
            assert G.N1==G.N2, "Const only applies to square G"
            C = C*numpy.identity(G.N1) 
        if C.shape !=(G.N1,G.N2) : raise RuntimeError, "Size of constant incorrect"

        G._tail.zero()
        G._tail[0].array[:,:] = C
        
        Function(lambda om : C, None)(G)
        return G
    
#########################################################################

class A_Omega_Plus_B(Base):
    def __init__ (self, A=1, B=0, Invert= False) :
        if isinstance(A,ArrayViewWithIndexConverter) : A=A.array[:,:] # trasnform to numpy.array
        if isinstance(B,ArrayViewWithIndexConverter) : B=B.array[:,:]# trasnform to numpy.array
        Base.__init__(self, A=A, B=B,Invert=Invert)
         
    def __call__(self,G) :
        A,B = self.A, self.B
        if G.mesh.TypeGF not in [GF_Type.Imaginary_Frequency, GF_Type.Real_Frequency] : 
            raise TypeError, "This initializer is only correct in frequency"

        if not isinstance(A,numpy.ndarray) : A = A*numpy.identity(G.N1) 
        if not isinstance(B,numpy.ndarray) : B = B*numpy.identity(G.N1) 
        if A.shape !=(G.N1,G.N2) : raise RuntimeError, "Size of A incorrect"
        if B.shape !=(G.N1,G.N2) : raise RuntimeError, "Size of B incorrect"

        t = G._tail
        t.zero()
        t[-1].array[:,:] = A
        t[0].array[:,:] = B

        Function(lambda om : A*om + B, None)(G)

        if self.Invert : G.invert()
        return G
#######################################

class OneFermionInTime(Base):
    def __init__ (self, l =0) :
         Base.__init__(self, L=l)
         
    def __call__(self,G) :
        L = self.L
        if G.mesh.TypeGF not in [GF_Type.Imaginary_Time] : 
            raise TypeError, "This initializer is only correct in frequency"

        t = G._tail
        t.zero()
        t[1].array[:,:] = 1
        t[2].array[:,:] = L
        t[3].array[:,:] = L*L
        
        fact = -1/(1+exp(-L*G.beta))
        Function(lambda t : fact* exp(-L*t), None)(G)
        return G


##################################################

def _SemiCircularDOS(half_bandwidth):
    """
       Semi_Circular DOS function
       Input : the 1/2 bandwidth
       Returns : a function omega-> dos(omega)
    """
    from math import sqrt,pi
    larg = half_bandwidth
    def semi(x):
        if (abs(x)<larg) : return sqrt( 1 - (x/larg)**2 )*2/pi/larg
        else: return 0.0
    return semi

def semi(x) :
    return _SemiCircularDOS(x)

##################################################

class SemiCircular (Base):
    def __init__ (self, half_bandwidth) :
        Base.__init__(self, half_bandwidth=half_bandwidth)

    def __call__(self,G) :
        D= self.half_bandwidth
        Id = numpy.identity(G.N1,numpy.complex_)
        if G.mesh.TypeGF == GF_Type.Imaginary_Frequency:
            f = lambda om : (om  - 1j*sign(om.imag)*sqrt(abs(om)**2 +  D*D))/D/D*2*Id
        elif G.mesh.TypeGF == GF_Type.Real_Frequency: 
            dos = _SemiCircularDOS (D)
            f = lambda om : -1j*pi*dos(om)*Id
        else :
            raise TypeError, "This initializer is only correct in frequency"

        t =G._tail
        t.zero()
        for i in range(G.N1):
            t[1].array[i,i] = 1.0
            t[3].array[i,i] = D**2/4.0
            t[5].array[i,i] = D**4/8.0

        #t[1].array[:,:] = 1.0
        #t[3].array[:,:] = D**2/4.0
        #t[5].array[:,:] = D**4/8.0

        t.changeOrderMax(5) # expansion is not valid above this
        
        #Mom[2,0,0] = 1.0
        #Mom[4,0,0] = D**2/4.0
        #Mom[6,0,0] = D**4/8.0

        Function(f,None)(G)
        return G

##################################################

class Wilson (Base):
    def __init__ (self, half_bandwidth) :
        Base.__init__(self, half_bandwidth=half_bandwidth)

    def __call__(self,G) :

        D = self.half_bandwidth
        Id = numpy.identity(G.N1,numpy.complex_)

        if G.mesh.TypeGF == GF_Type.Imaginary_Frequency:
            f = lambda om : (-1/(2.0*D)) * numpy.log((om-D)/(om+D))
        elif G.mesh.TypeGF == GF_Type.Real_Frequency: 
            def f(om):
              if (om > -D) and (om < D):
                return -numpy.log(abs(om-D)/abs(om+D))*Id/(2*D) - 1j*pi*Id/(2*D)
              else:
                return -numpy.log(abs(om-D)/abs(om+D))*Id/(2*D)
        else :
            raise TypeError, "This initializer is only correct in frequency"


        t = G._tail
        t.zero()
        for i in range(G.N1):
            t[1].array[i,i] = 1.0
            t[3].array[i,i] = D**2/3.0
            t[5].array[i,i] = D**4/5.0

        t.changeOrderMax(5) # expansion is not valid above this
        
        Function(f,None)(G)
        return G



"""
#########################################################################

class FreeBath (Function):

    ""@ Usage : 
         G <<= GFInitializer__to_FreeBath(eps0,t_array,eps_array)
       where : 
         * eps0 : blab la
         * t_array : 
         * eps_array : 
       Action : 
         Puts G to (omega - eps0 - sum_i t_array[i]^2/(omega - eps_array[i])
    ""@

    def __init__ (self,eps0,t_array,eps_array) :
        ""@    
        ""@
        Id = numpy.identity(G.Nnd,numpy.complex_)
        t,eps = [],[]
        dot = numpy.dot
	
        # check the types
        for t_ in t_array :
            if type(t_) != type(Id) :
                t.append(t_*Id)
            else :
                t.append(t_)
        for ep in eps_array  :
            if type(ep) != type(Id) :
                eps.append(ep*Id)
            else :
                eps.append(ep)
        #eps0= numpy.array(eps0)
	if type(eps0) != type(Id) :
            eps0= eps0*Id
        #
        def f(omega) : 
            res = 0*Id
            for tt,e in zip(t,eps) :
		res += dot(tt,dot(numpy.LinearAlgebra.inverse(omega*Id - e),numpy.transpose(tt)))
	    res = omega*Id - eps0 - res
	    return res

        # compute the moments
        tail = TailGf(0,4, []) ##### FAUX @@@@@@@@@@@@
        Mom = numpy.zeros((self.alpha_Max+2,self.Nnd,self.Nnd),numpy.complex_)
        Mom[0,:,:] = Id # 1/omega
        Mom[1,:,:] = - eps0
	Mom[2,:,:] = - sum([ dot(tt,tt) for tt in t])
        aux = [ dot(e,tt) for tt,e in zip(t,eps)]
        Mom[3,:,:] = - sum([ dot(tt,a) for tt,a in zip(t,aux)])
        Mom[4,:,:] = - sum([ dot(tt,dot(a,a)) for tt,a in zip(t,aux)])

        Function.__init__(self,f,Mom) # !!!!!!!!!! TRANSRCIP

#########################################################################


class FreeBath2 (Function):

    ""@ Usage : 
         G <<= GFInitializer__to_FreeBath2(eps0,t_array,eps_array)
       where : 
         * eps0 : blab la
         * t_array : 
         * eps_array : 
       Action : 
         Puts G to (omega - eps0 - sum_i t_array[i]^2/(omega - eps_array[i])
    ""@

    # --------------------------------------------------------------------

    def __init__ (self,eps0,t_array,eps_array) :
        ""@ ""@
        NBath = len(eps_array)

        Id_1 = numpy.identity(G.Nnd,numpy.complex_)
        Id_2 = numpy.identity(NBath,numpy.complex_)

        dot = numpy.dot
        if type(eps0) != type(Id_1) :
          eps0= eps0*Id_1

        e = 0*Id_2
        for i in xrange(NBath):
          e[i,i] = eps_array[i]

        def f(omega) : 
          res = dot(t_array,dot(numpy.LinearAlgebra.inverse(omega*Id_2 - e),numpy.transpose(t_array)))
          res = omega*Id_1 - eps0 - res
          return res

        # compute the moments
        Mom = numpy.zeros((self.alpha_Max+2,self.Nnd,self.Nnd),numpy.complex_)
        Mom[0,:,:] = Id_1
        Mom[1,:,:] = - eps0
        Mom[2,:,:] = - dot(t_array, numpy.transpose(t_array))
        Mom[3,:,:] = - dot(t_array, dot(e, numpy.transpose(t_array)))
        Mom[4,:,:] = - dot(t_array, dot(e, dot(e, numpy.transpose(t_array))))

        Function.__init__(self,f,Mom) # !!!!!!!!!! TRANSRCIP


#########################################################################

########################################################################


class SemiCircular_Nambu (Function):
    ""@ Usage : 
         G <<= GFInitializer__to_FreeBath2(D,h=0,mu=0,delta=0)
       where : 
         * D: half bandwidth 
         * h = 
         * mu = 
         * delta = 
       Action : 
       Sets G to a semi-circular density of states of width D

       TO DO : repair the h in matsubara !!!!
    ""@

    # --------------------------------------------------------------------

    def __init__ (self,D,h=0,mu=0,delta=0) :
	"
        if self.Nnd == 2:
          def f(om) :

            res = numpy.zeros((self.Nnd,self.Nnd),numpy.complex_)

            res[0,0] = (-1/D**2)*(-2*(h+mu+om)\
              +(1/numpy.sqrt((h+om)**2-delta**2))*(1j*(h+om-numpy.sqrt((h+om)**2-delta**2))*sign(om.imag)*numpy.sqrt(D**2-(mu-numpy.sqrt((h+om)**2-delta**2))**2))\
              +(1/numpy.sqrt((h+om)**2-delta**2))*(1j*(h+om+numpy.sqrt((h+om)**2-delta**2))*sign(om.imag)*numpy.sqrt(D**2-(mu+numpy.sqrt((h+om)**2-delta**2))**2)))

            res[1,1] = (-1/D**2)*(-2*(h-mu+om)\
              +(1/numpy.sqrt((h+om)**2-delta**2))*(1j*(h+om+numpy.sqrt((h+om)**2-delta**2))*sign(om.imag)*numpy.sqrt(D**2-(mu-numpy.sqrt((h+om)**2-delta**2))**2))\
              +(1/numpy.sqrt((h+om)**2-delta**2))*(1j*(h+om-numpy.sqrt((h+om)**2-delta**2))*sign(om.imag)*numpy.sqrt(D**2-(mu+numpy.sqrt((h+om)**2-delta**2))**2)))

            res[0,1] = (1/(D**2*numpy.sqrt((h+om)**2-delta**2)))*delta*\
              (-2*numpy.sqrt((h+om)**2-delta**2)+1j*sign(om.imag)*numpy.sqrt(D**2-(mu-numpy.sqrt((h+om)**2-delta**2))**2)\
                                              +1j*sign(om.imag)*numpy.sqrt(D**2-(mu+numpy.sqrt((h+om)**2-delta**2))**2))

            res[1,0] = res[0,1]

            return res

          Mom = numpy.zeros((self.alpha_Max+2,self.Nnd,self.Nnd),numpy.complex_)
          Mom[2,0,0] = 1.0
          Mom[3,0,0] = -h-mu
          Mom[4,0,0] = D**2/4.0+delta**2+(h+mu)**2
          Mom[5,0,0] = (1.0/4.0)*(-3*D**2*(h+mu)-4*((h+mu)**3+delta**2*(3*h+mu)))
          Mom[6,0,0] = D**4/8.0+delta**4+(h+mu)**4+2*delta**2*(3*h**2+2*h*mu+mu**2)+ \
                        (1.0/2.0)*D**2*(delta**2+3*(h+mu)**2)

          Mom[2,1,1] = 1.0
          Mom[3,1,1] = -h+mu
          Mom[4,1,1] = D**2/4.0+delta**2+(h-mu)**2
          Mom[5,1,1] = (1.0/4.0)*(-3*D**2*(h-mu)-4*((h-mu)**3+delta**2*(-3*h+mu)))
          Mom[6,1,1] = D**4/8.0+delta**4+(h-mu)**4+2*delta**2*(3*h**2-2*h*mu+mu**2)+ \
                        (1.0/2.0)*D**2*(delta**2+3*(h-mu)**2)

          Mom[2,0,1] = 0.0
          Mom[3,0,1] = -delta
          Mom[4,0,1] = 2.0*delta*h
          Mom[5,0,1] = -(1.0/4.0)*delta*(D**2+4*(delta**2+3*h**2+mu**2))
          Mom[6,0,1] = delta*h*(D**2+4*(delta**2+h**2+mu**2))

          Mom[2,1,0] = Mom[2,0,1]
          Mom[3,1,0] = Mom[3,0,1]
          Mom[4,1,0] = Mom[4,0,1]
          Mom[5,1,0] = Mom[5,0,1]
          Mom[6,1,0] = Mom[6,0,1]

          self.computeFrom(f,Mom)

    #-----------------------------------------------------
	  
"""

