
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

from pytriqs.applications.impurity_solvers.hubbard_I import Solver
import pytriqs.utility.mpi as mpi

class SolverOld(Solver):
    """
       Old construction of the Hubbard I Solver
    """

    def __init__(self, Beta, Uint, JHund, l, Nmsb=1025, T=None, UseSpinOrbit=False, Verbosity=0):


        Solver.__init__(self, beta=Beta, l=l, n_msb=Nmsb, use_spin_orbit=UseSpinOrbit)

        self.Uint = Uint
        self.JHund = JHund
        self.T = T
        self.Verbosity = Verbosity

        msg = """
**********************************************************************************
 Warning: You are using the old constructor for the solver. Beware that this will
 be deprecated in future versions. Please check the documentation.
**********************************************************************************
"""
        mpi.report(msg)


    def Solve(self, Iteration_Number=1, Test_Convergence=0.0001):

        self.solve(U_int = self.Uint, J_hund=self.JHund, T=self.T, verbosity=self.Verbosity,
                   Iteration_Number=Iteration_Number, Test_Convergence=Test_Convergence)
