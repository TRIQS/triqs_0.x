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

r""" """

import numpy
from math import *
from gf import MeshImFreq,TailGf #, MeshRealFrequency
from pytriqs.base.utility.myUtils import sign
from lazy_expressions import lazy_expr_terminal, transform, lazy_expr

def is_lazy(y) :
    #return type(y) in [ Omega_, lazy_expr]
    return  isinstance(y,(Omega_, lazy_expr,lazy_expr_terminal))

def is_scalar (x) : 
    return type(x) in [ type(1), type(1.0), type(1j), numpy.ndarray, numpy.int, numpy.int_, numpy.int8, numpy.int16, numpy.int32, numpy.float, numpy.float_, numpy.float32, numpy.float64, numpy.complex, numpy.complex_, numpy.complex64, numpy.complex128 ]

def convert_scalar_to_Const(expr) : 

  # if the expression is a pure scalar, replace it by Const
  t= expr.get_terminal()
  if is_scalar(t) : return lazy_expr( Const(t) )

  # otherwise : replace all scalar appearing in +/- operations by Const
  def act (tag, childs) : 
        if tag in ["+", "-"] :
            for n,c in enumerate(childs) : 
                t = c.get_terminal()
                if is_scalar(t) : childs[n] =  Const (t)
        return (tag,childs)

  return transform(expr, act)
    
#########################################################################

class base (lazy_expr_terminal) :
    def __init__(self,**kargs) :
        self.__dict__.update(kargs)
                
#########################################################################

class Function (base):
    r"""   
       Stores a python function and a tail.
       
       If the Green's function is defined on an array of points :math:`x_i`, then it will be initialized to :math:`F(x_i)`.
    """
    def __init__ (self, F, Tail=None) : 
        r"""
          :param F: the function  :math:`\omega \rightarrow F(\omega)`
          :param Tail: The tail. Use None if you don't use any tail (will be put to 0)
        """
        base.__init__(self,F=F, Tail=Tail)
        
    def __call__(self,G) :
        if not(callable(self.F)) : raise RuntimeError, "GFInitializer.Function : f must be callable"
        res = G.data[:,:,:]
        try :
            for n,om in enumerate(G.mesh) : res[:,:,n] = self.F(om)
        except :
            raise RuntimeError, "The given function has a problem..."
        if self.Tail : G.tail.copy_from(self.Tail)
        return G

#########################################################################

class Const(base):
    def __init__ (self, C) :
        base.__init__(self, C=C)
         
    def __call__(self,G) :
        C = self.C
        if G.mesh.__class__.__name__ not in ['MeshImFreq', 'MeshRealFrequency']:
            raise TypeError, "This initializer is only correct in frequency"

        if not isinstance(C,numpy.ndarray) : 
            assert G.N1==G.N2, "Const only applies to square G"
            C = C*numpy.identity(G.N1) 
        if C.shape !=(G.N1,G.N2) : raise RuntimeError, "Size of constant incorrect"

        t = G.tail
        G.tail = TailGf(shape = t.shape, size = t.size, order_min=0)
        G.tail[0][:,:] = C
        
        Function(lambda om : C, None)(G)
        return G
    
#########################################################################

class Omega_(base):
    r"""The function :math:`\omega \rightarrow \omega` """
    def __str__(self) : return "Omega" 
    def __call__(self,G) :
        if G.mesh.__class__.__name__ not in ['MeshImFreq', 'MeshRealFrequency']:
            raise TypeError, "This initializer is only correct in frequency"

        t,Id  = G.tail, numpy.identity(G.N1)
        G.tail = TailGf(shape = t.shape, size = t.size, order_min=-1)
        G.tail[-1][:,:] = Id
        Function(lambda om : om*numpy.identity(G.N1), None)(G)
        return G

Omega = Omega_()
iOmega_n = Omega_()

##########################################################################

class A_Omega_Plus_B(base):
    "deprecated. do not use"
    def __init__ (self, A=1, B=0, Invert= False) :
        base.__init__(self, A=A, B=B,Invert=Invert)
         
    def __call__(self,G) :
        A,B = self.A, self.B
        if G.mesh.__class__.__name__ not in ['MeshImFreq', 'MeshRealFrequency']:
            raise TypeError, "This initializer is only correct in frequency"

        if not isinstance(A,numpy.ndarray) : A = A*numpy.identity(G.N1) 
        if not isinstance(B,numpy.ndarray) : B = B*numpy.identity(G.N1) 
        if A.shape !=(G.N1,G.N2) : raise RuntimeError, "Size of A incorrect"
        if B.shape !=(G.N1,G.N2) : raise RuntimeError, "Size of B incorrect"

        t,Id  = G.tail, numpy.identity(G.N1)
        G.tail = TailGf(shape = t.shape, size = t.size, order_min=-1)
        G.tail[-1][:,:] = A
        G.tail[0][:,:] = B
        
        Function(lambda om : A*om + B, None)(G)

        if self.Invert : G.invert()
        return G

#######################################

class OneFermionInTime(base):
    def __init__ (self, l =0) :
         base.__init__(self, L=l)
         
    def __call__(self,G) :
        L = self.L
        if G.mesh.TypeGF not in [GF_Type.Imaginary_Time] : 
            raise TypeError, "This initializer is only correct in frequency"

        t,Id  = G.tail, numpy.identity(G.N1)
        G.tail = TailGf(shape = t.shape, size = t.size, order_min=1)
        t[1][:,:] = 1
        t[2][:,:] = L
        t[3][:,:] = L*L
        
        fact = -1/(1+exp(-L*G.Beta))
        Function(lambda t : fact* exp(-L*t) *Id, None)(G)
        return G


##################################################

def _SemiCircularDOS(HalfBandwidth):
    """
       Semi_Circular DOS function
       Input : the 1/2 bandwidth
       Returns : a function omega-> dos(omega)
    """
    from math import sqrt,pi
    larg = HalfBandwidth
    def semi(x):
        if (abs(x)<larg) : return sqrt( 1 - (x/larg)**2 )*2/pi/larg
        else: return 0.0
    return semi

def semi(x) :
    return _SemiCircularDOS(x)

##################################################

class SemiCircular (base):
    r"""Hilbert transform of a semi circular density of state, i.e.

     .. math ::
        g(z) = \int \frac{A(\omega)}{z-\omega} d\omega
        
    where :math:`A(\omega) = \theta( D - |\omega|) 2 \sqrt{ D^2 - \omega^2}/(\pi D^2)` 
      
     (only works in combination with frequency Green's functions).
    """
    def __init__ (self, HalfBandwidth) :
        """ :param HalfBandwidth:  :math:`D`, the half bandwidth of the semicircular"""
        base.__init__(self, HalfBandwidth=HalfBandwidth)

    def __str__(self) : return "SemiCircular(%s)"%self.HalfBandwidth 

    def __call__(self,G) :
        D= self.HalfBandwidth
        Id = numpy.identity(G.N1,numpy.complex_)
        if type(G.mesh) == MeshImFreq:
            f = lambda om : (om  - 1j*sign(om.imag)*sqrt(abs(om)**2 +  D*D))/D/D*2*Id
        elif type(G.mesh) == MeshRealFrequency:
            def f(om):
              if (om > -D) and (om < D):
                return (2.0/D**2) * (om - 1j* sqrt(D**2 - om**2))
              else:
                return (2.0/D**2) * (om - sign(om) * sqrt(om**2 - D**2))
        else :
            raise TypeError, "This initializer is only correct in frequency"

        # Let's create a new tail
        t = G.tail
        G.tail = TailGf(shape = t.shape, size = t.size, order_min=1)
        for i in range(G.N1):
            G.tail[1][i,i] = 1.0
            G.tail[3][i,i] = D**2/4.0
            G.tail[5][i,i] = D**4/8.0
 
        Function(f,None)(G)
        return G

##################################################

class Wilson (base):
    r"""The Hilbert transform of a flat density of states, with cut-off

    .. math ::
        g(z) = \int \frac{A(\omega)}{z-\omega} d\omega
        
    where :math:`A(\omega) = \theta( D^2 - \omega^2)/(2D)` 
      
    (only works in combination with frequency Green's functions).
    """
    def __init__ (self, HalfBandwidth) :
        """:param HalfBandwidth: :math:`D`, the half bandwidth """
        base.__init__(self, HalfBandwidth=HalfBandwidth)

    def __str__(self) : return "Wilson(%s)"%HalfBandwidth 

    def __call__(self,G) :

        D = self.HalfBandwidth
        Id = numpy.identity(G.N1,numpy.complex_)

        if type(G.mesh) == MeshImFreq:
            f = lambda om : (-1/(2.0*D)) * numpy.log((om-D)/(om+D)) * Id
        elif type(G.mesh) == MeshRealFrequency:
            def f(om):
              if (om > -D) and (om < D):
                return -numpy.log(abs(om-D)/abs(om+D))*Id/(2*D) - 1j*pi*Id/(2*D)
              else:
                return -numpy.log(abs(om-D)/abs(om+D))*Id/(2*D)
        else :
            raise TypeError, "This initializer is only correct in frequency"

        # Let's create a new tail
        G.tail = TailGf(shape = G.tail.shape, size=5, order_min=1)
        for i in range(G.N1):
            G.tail[1][i,i] = 1.0
            G.tail[3][i,i] = D**2/3.0
            G.tail[5][i,i] = D**4/5.0

        Function(f,None)(G)
        return G


##################################################

class Fourier (base):
    r"""
    The Fourier transform as a lazy expression
    """
    def __init__ (self, G) :
        """:param G: :math:`G`, the function to be transformed. Must in the time domain"""
        base.__init__(self, G = G)

    def __str__(self) : return "Fourier(%s)"%self.G.Name

    def __call__(self,G2) :
        G2.set_from_fourier_of(self.G)
        return G2

class InverseFourier (base):
    r"""
    The Inverse Fourier transform as a lazy expression
    """
    def __init__ (self, G) :
        """:param G: :math:`G`, the function to be transformed. Must in the frequency domain"""
        base.__init__(self, G = G)

    def __str__(self) : return "InverseFourier(%s)"%self.G.Name

    def __call__(self,G2) :
        G2.set_from_inverse_fourier_of(self.G)
        return G2

class LegendreToMatsubara (base):
    r"""
    The transformation from Legendre to Matsubara as a lazy expression
    """
    def __init__ (self, G) :
        """:param G: :math:`G`, the function to be transformed. Must in the Legendre domain"""
        base.__init__(self, G = G)

    def __str__(self) : return "LegendreToMatsubara(%s)"%self.G.Name

    def __call__(self,G2) :
        G2.setFromLegendre(self.G)
        return G2

class MatsubaraToLegendre (base):
    r"""
    The transformation from Legendre to Matsubara as a lazy expression
    """
    def __init__ (self, G) :
        """:param G: :math:`G`, the function to be transformed. Must in the Matsubara domain"""
        base.__init__(self, G = G)

    def __str__(self) : return "MatsubaraToLegendre(%s)"%self.G.Name

    def __call__(self,G2) :
        G2.setFromMatsubara(self.G)
        return G2

