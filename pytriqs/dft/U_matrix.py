
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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

# calculates the four index U matrix

import numpy
from types import *
from math import sqrt
import copy
from pytriqs.dft.vertex import u4ind

class Umatrix:
    """calculates, stores, and manipulates the four index U matrix"""

    def __init__(self, l, U_interact=0, J_Hund=0):

        self.l = l
        self.U_av = U_interact
        self.J = J_Hund

        self.N = 2*l+1         # multiplicity

        #self.Ucmplx = numpy.zeros([self.N,self.N,self.N,self.N],numpy.float_)
        #self.Ucubic = numpy.zeros([self.N,self.N,self.N,self.N],numpy.float_)


    def __call__(self, T = None, Rcl = None):
        """calculates the four index matrix. Slater parameters can be provided in Rcl, 
           and a transformation matrix from complex harmonics to a specified other representation (e.g. cubic). 
           If T is not given, use standard complex harmonics."""

        if Rcl is None: Rcl = self.getRcl(self.U_av,self.J,self.l)

        if (T is None): 
            TM = numpy.identity(self.N,numpy.complex_)
        else:
            TM = T
        self.Nmat = len(TM)

        self.Ufull = u4ind(Rcl,TM)
        
                    
        
    def ReduceMatrix(self):
        """Reduces the four-index matrix to two-index matrices."""

        if (self.N==self.Nmat):
            self.U  = numpy.zeros([self.N,self.N],numpy.float_)    # matrix for same spin
            self.Up = numpy.zeros([self.N,self.N],numpy.float_)    # matrix for opposite spin

            for m in range(self.N):
                for mp in range(self.N):
                    self.U[m,mp]  = self.Ufull[m,mp,m,mp].real - self.Ufull[m,mp,mp,m].real
                    self.Up[m,mp] = self.Ufull[m,mp,m,mp].real
        else:
            self.U  = numpy.zeros([self.Nmat,self.Nmat],numpy.float_)    # matrix 
            for m in range(self.Nmat):
                for mp in range(self.Nmat):
                    self.U[m,mp]  = self.Ufull[m,mp,m,mp].real - self.Ufull[m,mp,mp,m].real
           
        
           


    def getRcl(self,Uint,JH,l):

        #Rcl = numpy.array([0.0, 0.0, 0.0, 0.0],numpy.float_)
        xx = l+1
        Rcl = numpy.zeros([xx],numpy.float_)
        if(l==2):
            Rcl[0] = Uint
            Rcl[1] = JH * 14.0 / (1.0 + 0.63)
            Rcl[2] = 0.630 * Rcl[1]
        elif(l==3):
            Rcl[0] = Uint
            Rcl[1] = 6435.0 * JH / (286.0 + 195.0 * 451.0 / 675.0 + 250.0 * 1001.0 / 2025.0)
            Rcl[2] = 451.0 * Rcl[1] / 675.0
            Rcl[3] = 1001.0 * Rcl[1] / 2025.0

        return Rcl
