
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

from math import *
import numpy
from pytriqs.base.gf_local import BlockGf
from pytriqs.base.utility.my_utils import sum_list
from pytriqs.solvers.operators import *
from pytriqs.solvers.ctqmc_hyb import Solver

#########################################
#
#  Solver for 2x2 plaquette ROTATION
#
#########################################

class Solver_2x2_Para_Hubbard (Solver) : 
   def __init__(self,Beta,U_interact, **kw) : 
      
      self.P = numpy.array([[1,1,-1,-1], [1,-1,1j,-1j],[1,-1,-1j,1j],[1,1,1,1]])/2
      self.Pinv = numpy.linalg.inv(self.P)

      def Cdiag(s):
         s= '-%s'%s
         return [ C('+1'+s,1), C('-1'+s,1), C('+i'+s,1),C('-i'+s,1) ] 

      Cnat,Nnat={},{}
      for i in range(4) :
         s= 'up-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('up'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]
         s= 'down-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('down'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]

      def symm(s):
         s= '-%s'%s
         Cdag('+1'+s,1).symmetry()['Z4'] = 1
         Cdag('-1'+s,1).symmetry()['Z4'] = -1
         Cdag('+i'+s,1).symmetry()['Z4'] = 1j
         Cdag('-i'+s,1).symmetry()['Z4'] = -1j
      symm('up'); symm('down')

      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact* sum_list( [ Nnat['up-%d'%(i+1)]* Nnat['down-%d'%(i+1)]  for i in range (4)]) #!!

      N_u = sum_list( [ Nnat['up-%d'%(i+1)] for i in range(4) ] )
      N_d = sum_list( [ Nnat['down-%d'%(i+1)] for i in range(4) ] )

      Quantum_Numbers = { 'Z4' : 'Z4', 'N_u' : N_u, 'N_d' : N_d }

      Solver.__init__(self,Beta=Beta,GFstruct=[ ('+1-up',[1]), ('-1-up',[1]), ('+i-up',[1]), ('-i-up',[1]),
                                           ('+1-down',[1]), ('-1-down',[1]), ('+i-down',[1]), ('-i-down',[1]) ], 
                      H_Local = Hamiltonian,Quantum_Numbers= Quantum_Numbers )
      self.N_Cycles  = 10

   def Transform_RealSpace_to_SymmetryBasis(self,IN, OUT = None):
      """ IN[i,j] in real space --> OUT into the symmetry basis"""
      OUT = OUT if OUT else BlockGf(G)
      for sig,B in IN:
         for k,ind in enumerate(['+1-', '-1-','+i-', '-i-']) : 
            OUT[ind+sig] = sum_list ( [ sum_list ( [ B[i+1,j+1]* self.P[j,k] * self.Pinv[k,i] for j in range(4) ] ) for i in range(4) ] ) 
      return OUT

   def Transform_SymmetryBasis_toRealSpace(self,IN, OUT = None):
      """  IN : in symmetry cluster indices. Returns OUT to real space"""
      OUT = OUT if OUT else BlockGf(G)
      for sig,B in OUT : 
         for i in range(4):
            for j in range(4):
               B[i+1,j+1]  = sum_list( [ self.P[i,k]* self.Pinv[k,j]* IN[ind+sig] 
                                    for k,ind in enumerate(['+1-', '-1-','+i-', '-i-']) ])
      return OUT

#########################################
#
#  Solver for 2x2 plaquette MOMENTUM
#
#########################################

class Solver_2x2_Para_Hubbard_Momentum (Solver) : 
   def __init__(self,Beta,U_interact) : 
      
      self.P = numpy.array([[1,1,1,1], [1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])/2.
      self.Pinv = numpy.linalg.inv(self.P)

      def Cdiag(s):
         s= '-%s'%s
         return [ C('00'+s,1), C('10'+s,1), C('01'+s,1),C('11'+s,1) ] 

      Cnat,Nnat={},{}
      for i in range(4) :
         s= 'up-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('up'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]
         s= 'down-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('down'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]

      def symm(s):
         s= '-%s'%s
         Cdag('00'+s,1).symmetry()['kx'] =  1
         Cdag('10'+s,1).symmetry()['kx'] = -1
         Cdag('01'+s,1).symmetry()['kx'] =  1
         Cdag('11'+s,1).symmetry()['kx'] = -1
         Cdag('00'+s,1).symmetry()['ky'] =  1
         Cdag('10'+s,1).symmetry()['ky'] =  1
         Cdag('01'+s,1).symmetry()['ky'] = -1
         Cdag('11'+s,1).symmetry()['ky'] = -1
      symm('up'); symm('down')

      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact* sum_list( [ Nnat['up-%d'%(i+1)]* Nnat['down-%d'%(i+1)]  for i in range (4)]) #!!

      N_u = sum_list( [ Nnat['up-%d'%(i+1)] for i in range(4) ] )
      N_d = sum_list( [ Nnat['down-%d'%(i+1)] for i in range(4) ] )

      Quantum_Numbers = { 'kx' : 'kx', 'ky' : 'ky', 'N_u' : N_u, 'N_d' : N_d }

      Solver.__init__(self,Beta=Beta,GFstruct=[ ('00-up',[1]), ('10-up',[1]), ('01-up',[1]), ('11-up',[1]),
                                                ('00-down',[1]), ('10-down',[1]), ('01-down',[1]), ('11-down',[1]) ],
                      H_Local = Hamiltonian,Quantum_Numbers= Quantum_Numbers )

      self.N_Cycles  = 10

##############################
#
#  Solver for a 3-band model
#
###############################

class Solver_3_band_Hubbard (Solver) : 
   def __init__(self,Beta,U_interact,J_Hund) : 
     
      GFstruct = [ ('1up',[1]), ('1down',[1]), ('23up',[1,2]), ('23down',[1,2]) ]
	 
      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      #Hamiltonian =   U_interact * sum_list( [ N('%d-up'%(i+1),1)*N('%d-down'%(i+1),1)  for i in range(3) ] )
      Hamiltonian =   U_interact * ( N('1up',1)*N('1down',1) + N('23up',1)*N('23down',1) + N('23up',2)*N('23down',2) )
      Hamiltonian +=  (U_interact-2.0*J_Hund) * ( N('1up',1) * (N('23down',1)+N('23down',2)) + N('23up',1) * N('23down',2) )
      Hamiltonian +=  (U_interact-2.0*J_Hund) * ( N('1down',1) * (N('23up',1)+N('23up',2)) + N('23up',2) * N('23down',1) )
      for s in ['up','down']:
        Hamiltonian += (U_interact-3.0*J_Hund) * (N('1%s'%s,1) * (N('23%s'%s,1)+N('23%s'%s,2)) + N('23%s'%s,1)*N('23%s'%s,2) )

#      for i in range(2):
#        for j in range(i+1,3):
#          Hamiltonian += (U_interact-2*J_Hund) * ( N('%d-up'%(i+1),1) * N('%d-down'%(j+1),1) +
#                                                   N('%d-down'%(i+1),1) * N('%d-up'%(j+1),1) ) + \
#                         (U_interact-3*J_Hund) * ( N('%d-up'%(i+1),1) * N('%d-up'%(j+1),1) +
#                                                   N('%d-down'%(i+1),1) * N('%d-down'%(j+1),1) )

#      Quantum_Numbers = { 'N1up' : N('1-up',1), 'N1down' : N('1-down',1),
#                         'N2up' : N('2-up',1), 'N2down' : N('2-down',1),
#                         'N3up' : N('3-up',1), 'N3down' : N('3-down',1) }
      Quantum_Numbers = { 'N1up' : N('1up',1), 'N1down' : N('1down',1),
                         'N23up' : N('23up',1)+N('23up',2),
                         'N23down' : N('23down',1)+N('23down',2) }

      Solver.__init__(self,
                      Beta = Beta,
                      GFstruct = GFstruct,
                      H_Local = Hamiltonian,
                      Quantum_Numbers = Quantum_Numbers )
      self.N_Cycles  = 10000

##############################
#
#  Solver for a 2- or 3-band model
#
###############################

class Solver_2_or_3_bands_Hubbard (Solver) : 
   def __init__(self,Beta,Norb,U_interact,J_Hund) : 
      
      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact * sum_list( [ N('up',i)*N('down',i) for i in range(Norb) ] )
      for i in range(Norb-1):
        for j in range(i+1,Norb):
          Hamiltonian += (U_interact-2*J_Hund) * ( N('up',i) * N('down',j) +
                                                   N('down',i) * N('up',j) ) + \
                         (U_interact-3*J_Hund) * ( N('up',i) * N('up',j) +
                                                   N('down',i) * N('down',j) )

      #Quantum_Numbers = {}
      #for i in range(Norb):
      #  Quantum_Numbers['N%sup'%(i+1)]   = N('up',i)
      #  Quantum_Numbers['N%sdown'%(i+1)] = N('down',i)

      Ntot = sum_list( [ N(s,i) for s in ['up','down'] for i in range(Norb) ] )
      Sz = sum_list( [ N('up',i)-N('down',i) for i in range(Norb) ] )

      Quantum_Numbers = {'Ntot' : Ntot, 'Sz' : Sz}

      Solver.__init__(self,
                      Beta = Beta,
                      GFstruct = [ ('%s'%(ud),[n for n in range(Norb)]) for ud in ['up','down'] ],
                      H_Local = Hamiltonian,
                      Quantum_Numbers = Quantum_Numbers )
      self.N_Cycles  = 10000

##############################
#
#  Solver for the crystal field problem
#
###############################

class Solver_XField (Solver) : 
   def __init__(self,Beta,U_interact) : 


     Hamiltonian = U_interact*( N('1-up',1)*N('1-down',1) + N('1-up',1)*N('2-up',1) +
                                N('1-up',1)*N('2-down',1) + N('1-down',1)*N('2-up',1) +
                                N('1-down',1)*N('2-down',1) + N('2-up',1)*N('2-down',1) )

     Quantum_Numbers = { 'Nup'   : N('1-up',1)+N('2-up',1),
                        'Ndown' : N('1-down',1)+N('2-down',1),
                        'Lz'    : N('1-up',1)+N('1-down',1)-N('2-up',1)-N('2-down',1) }

     Solver.__init__(self,
                     Beta = Beta,
                     GFstruct = [ ('1-up',[1]), ('1-down',[1]), ('2-up',[1]), ('2-down',[1]) ],
                     H_Local = Hamiltonian,
                     Quantum_Numbers = Quantum_Numbers )
     self.N_Cycles  = 1000000
     self.Nmax_Matrix = 300

#########################################
#
#  Solver for the dimer
#
#########################################

class Solver_Dimer (Solver) : 
   def __init__(self,Beta,U_interact) : 
      
      C_up_1 = (C('Up+',1) + C('Up-',1))/sqrt(2)
      C_up_2 = (C('Up+',1) - C('Up-',1))/sqrt(2)
      C_do_1 = (C('Do+',1) + C('Do-',1))/sqrt(2)
      C_do_2 = (C('Do+',1) - C('Do-',1))/sqrt(2)

      Cdag('Up+',1).symmetry()['parity'] = 1
      Cdag('Up-',1).symmetry()['parity'] = -1
      Cdag('Do+',1).symmetry()['parity'] = 1
      Cdag('Do-',1).symmetry()['parity'] = -1

      N_up_1 = C_up_1.dagger()*C_up_1
      N_up_2 = C_up_2.dagger()*C_up_2
      N_do_1 = C_do_1.dagger()*C_do_1
      N_do_2 = C_do_2.dagger()*C_do_2

      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact * ( N_up_1*N_do_1  + N_up_2*N_do_2 )

      N_u = N('Up-',1)+N('Up+',1)
      N_d = N('Do-',1)+N('Do+',1)

      Quantum_Numbers = { 'parity' : 'parity', 'N_u' : N_u, 'N_d' : N_d }

      Solver.__init__(self,
                      Beta = Beta,
                      GFstruct = [ ('Up+',[1]), ('Up-',[1]), ('Do+',[1]), ('Do-',[1]) ],
                      H_Local = Hamiltonian,
                      Quantum_Numbers = Quantum_Numbers)

      self.N_Cycles  = 100000
      self.Length_Cycle = 100
      self.N_Frequencies_Accumulated = 120
      self.fitting_Frequency_Start = 100
      self.Use_Segment_Picture = False

#########################################
#
#  Solver for the Anderson model
#
#########################################

class Solver_Anderson (Solver) : 
   def __init__(self,Beta,U_interact) : 
      
      C_up_1 = (C('Up+',1) + C('Up-',1))/sqrt(2)
      C_up_2 = (C('Up+',1) - C('Up-',1))/sqrt(2)
      C_do_1 = (C('Do+',1) + C('Do-',1))/sqrt(2)
      C_do_2 = (C('Do+',1) - C('Do-',1))/sqrt(2)

      Cdag('Up+',1).symmetry()['parity'] = 1
      Cdag('Up-',1).symmetry()['parity'] = -1
      Cdag('Do+',1).symmetry()['parity'] = 1
      Cdag('Do-',1).symmetry()['parity'] = -1

      N_up_1 = C_up_1.dagger()*C_up_1
      N_up_2 = C_up_2.dagger()*C_up_2
      N_do_1 = C_do_1.dagger()*C_do_1
      N_do_2 = C_do_2.dagger()*C_do_2

      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact * ( N_up_1*N_do_1  + N_up_2*N_do_2 )

      N_u = N('Up-',1)+N('Up+',1)
      N_d = N('Do-',1)+N('Do+',1)

      Quantum_Numbers = { 'parity' : 'parity', 'N_u' : N_u, 'N_d' : N_d }

      Solver.__init__(self,
                      Beta = Beta,
                      GFstruct = [ ('Up+',[1]), ('Up-',[1]), ('Do+',[1]), ('Do-',[1]) ],
                      H_Local = Hamiltonian,
                      Quantum_Numbers = Quantum_Numbers,
                      Nmax = 2000)

      self.N_Cycles  = 100000
      self.N_Frequencies_Accumulated = 4*int(0.075*Beta/(2*3.1415))
      self.fitting_Frequency_Start = 3*int(0.075*Beta/(2*3.1415))
      self.Length_Cycle = 100
      self.Nmax_Matrix = 200
      self.Use_Segment_Picture = False

################################################
#
#  Solver for the Anderson model on a plaquette
#
################################################

class Solver_Anderson_4sites (Solver) : 
   def __init__(self,Beta,U_interact) : 

      self.P = numpy.array([[1,1,-1,-1], [1,-1,1j,-1j],[1,-1,-1j,1j],[1,1,1,1]])/2
      self.Pinv = numpy.linalg.inv(self.P)

      def Cdiag(s):
         s= '-%s'%s
         return [ C('+1'+s,1), C('-1'+s,1), C('+i'+s,1),C('-i'+s,1) ] 

      Cnat,Nnat={},{}
      for i in range(4) :
         s= 'up-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('up'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]
         s= 'down-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('down'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]

      def symm(s):
         s= '-%s'%s
         Cdag('+1'+s,1).symmetry()['Z4'] = 1
         Cdag('-1'+s,1).symmetry()['Z4'] = -1
         Cdag('+i'+s,1).symmetry()['Z4'] = 1j
         Cdag('-i'+s,1).symmetry()['Z4'] = -1j
      symm('up'); symm('down')

      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact* sum_list( [ Nnat['up-%d'%(i+1)]* Nnat['down-%d'%(i+1)] for i in range (4)])

      N_u = sum_list( [ Nnat['up-%d'%(i+1)] for i in range(4) ] )
      N_d = sum_list( [ Nnat['down-%d'%(i+1)] for i in range(4) ] )

      Quantum_Numbers = { 'Z4' : 'Z4', 'N_u' : N_u, 'N_d' : N_d }

      Solver.__init__(self,
                      Beta = Beta,
                      GFstruct = [ ('+1-up',[1]), ('-1-up',[1]), ('+i-up',[1]), ('-i-up',[1]),
                                   ('+1-down',[1]), ('-1-down',[1]), ('+i-down',[1]), ('-i-down',[1]) ], 
                      H_Local = Hamiltonian,
                      Quantum_Numbers = Quantum_Numbers,
                      Nmax = 2000)

      self.N_Cycles  = 100000
      self.N_Frequencies_Accumulated = 4*int(0.075*Beta/(2*3.1415))
      self.fitting_Frequency_Start = 3*int(0.075*Beta/(2*3.1415))
      self.Length_Cycle = 100
      self.Nmax_Matrix = 200
      self.Use_Segment_Picture = False

   def Transform_RealSpace_to_SymmetryBasis(self,IN, OUT = None):
      """ IN[i,j] in real space --> OUT into the symmetry basis"""
      OUT = OUT if OUT else BlockGf(G)
      for sig,B in IN:
         for k,ind in enumerate(['+1-', '-1-','+i-', '-i-']) : 
            OUT[ind+sig] = sum_list ( [ sum_list ( [ B[(i+1,1),(j+1,1)]* self.P[j,k] * self.Pinv[k,i] for j in range(4) ] ) for i in range(4) ] ) 
      return OUT

   def Transform_SymmetryBasis_to_RealSpace(self,IN, OUT = None):
      """  IN : in symmetry cluster indices. Returns OUT to real space"""
      OUT = OUT if OUT else BlockGf(G)
      for sig,B in OUT : 
         for i in range(4):
            for j in range(4):
               B[(i+1,1),(j+1,1)]  = sum_list( [ self.P[i,k]* self.Pinv[k,j]* IN[ind+sig] 
                                    for k,ind in enumerate(['+1-', '-1-','+i-', '-i-']) ])
      return OUT

################################################
#
#  Solver for the Anderson model on a plaquette MOMENTUM
#
################################################

class Solver_Anderson_4sites_Momentum (Solver) : 
   def __init__(self,Beta,U_interact) : 

      self.P = numpy.array([[1,1,1,1], [1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])/2.
      self.Pinv = numpy.linalg.inv(self.P)

      def Cdiag(s):
         s= '-%s'%s
         return [ C('00'+s,1), C('10'+s,1), C('01'+s,1),C('11'+s,1) ] 

      Cnat,Nnat={},{}
      for i in range(4) :
         s= 'up-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('up'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]
         s= 'down-%d'%(i+1)
         Cnat[s] = sum_list( [ self.P[i,j]*c  for j,c in enumerate(Cdiag('down'))] )
         Nnat[s] = Cnat[s].dagger() * Cnat[s]

      def symm(s):
         s= '-%s'%s
         Cdag('00'+s,1).symmetry()['kx'] =  1
         Cdag('10'+s,1).symmetry()['kx'] = -1
         Cdag('01'+s,1).symmetry()['kx'] =  1
         Cdag('11'+s,1).symmetry()['kx'] = -1
         Cdag('00'+s,1).symmetry()['ky'] =  1
         Cdag('10'+s,1).symmetry()['ky'] =  1
         Cdag('01'+s,1).symmetry()['ky'] = -1
         Cdag('11'+s,1).symmetry()['ky'] = -1
      symm('up'); symm('down')

      # NB : the Hamiltonian should NOT contain the quadratic part which is in G0
      Hamiltonian = U_interact* sum_list( [ Nnat['up-%d'%(i+1)]* Nnat['down-%d'%(i+1)] for i in range (4)])

      N_u = sum_list( [ Nnat['up-%d'%(i+1)] for i in range(4) ] )
      N_d = sum_list( [ Nnat['down-%d'%(i+1)] for i in range(4) ] )

      Quantum_Numbers = { 'kx' : 'kx', 'ky' : 'ky', 'N_u' : N_u, 'N_d' : N_d }

      Solver.__init__(self,
                      Beta = Beta,
                      GFstruct = [ ('00-up',[1]), ('10-up',[1]), ('01-up',[1]), ('11-up',[1]),
                                   ('00-down',[1]), ('10-down',[1]), ('01-down',[1]), ('11-down',[1]) ], 
                      H_Local = Hamiltonian,
                      Quantum_Numbers = Quantum_Numbers,
                      Nmax = 2000)

      self.N_Cycles  = 100000
      self.N_Frequencies_Accumulated = 4*int(0.075*Beta/(2*3.1415))
      self.fitting_Frequency_Start = 3*int(0.075*Beta/(2*3.1415))
      self.Length_Cycle = 100
      self.Nmax_Matrix = 200
      self.Use_Segment_Picture = False

   def Transform_RealSpace_to_SymmetryBasis(self,IN, OUT = None):
      """ IN[i,j] in real space --> OUT into the symmetry basis"""
      OUT = OUT if OUT else BlockGf(G)
      for sig,B in IN:
         for k,ind in enumerate(['00-', '10-','01-', '11-']) : 
            OUT[ind+sig] = sum_list ( [ sum_list ( [ B[(i+1,1),(j+1,1)]* self.P[j,k] * self.Pinv[k,i] for j in range(4) ] ) for i in range(4) ] ) 
      return OUT

   def Transform_SymmetryBasis_to_RealSpace(self,IN, OUT = None):
      """  IN : in symmetry cluster indices. Returns OUT to real space"""
      OUT = OUT if OUT else BlockGf(G)
      for sig,B in OUT : 
         for i in range(4):
            for j in range(4):
               B[(i+1,1),(j+1,1)]  = sum_list( [ self.P[i,k]* self.Pinv[k,j]* IN[ind+sig] 
                                    for k,ind in enumerate(['00-', '10-','01-', '11-']) ])
      return OUT
