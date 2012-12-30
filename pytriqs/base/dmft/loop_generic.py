
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

__all__ = ["DMFT_Loop_Generic"]

import shelve,os,signal,types
from pytriqs.solvers import *
import pytriqs.base.utility.MPI as MPI

class DMFT_Loop_Generic:
  """ This class provides a generic and simple DMFT loop for convenience.
      Given a list of solvers, its method `run` runs the generic DMFT loop.

      Usage : 

        * Derive and redefine Self_Consistency() (and maybe PostSolver)

      Useful elements : 

        * self.Iteration_Number is the number of the DMFT loop iteration
              
      """

  def __init__(self, Solver_List ,  **AdditionnalParameters) : 
    """

    :param Solver_List:  * List of the solvers to be called in the loop. 
                         *  If there is only one solver S, the syntax Solver_List = S 
                            is tolerated instead of Solver_List = [S]

    :param AdditionnalParameters: 
                        * all keywords arguments are simply put in the namespace
                        * of the loop e.g. :: 
                            
                             myloop(Solver_List = S,Chemical_potential= 2.0, theta = 78)
                          
                          myloop will be constructed with ::

                            self.Chemical_potential = 2.0
                            self.theta = 78

    """
    import warnings
    warnings.warn("The class %s is deprecated and scheduled for removal in later version of TRIQS. Please upgrade your code (write your own loop)"%(self.__class__.__name__), DeprecationWarning)
    self.SolverList = Solver_List if type(Solver_List) in [types.ListType,types.TupleType] else [Solver_List]
    for s in self.SolverList : assert isinstance(s,SolverBase), "Oops, one element of Solver_List is not a solver !!!"

    self.Iteration_Number  = 1
    # Put the additionnal parameters in my name space.
    self.__dict__.update(AdditionnalParameters)
     
  ##--------------------------------------------
      
  def PostSolver(self) : 
      """ Do nothing by default """
      pass
     
  ##-------------------------------------------

  # Obsolete : due to old RUPC config.
  def soft_kill_signal(self) : return os.path.exists('.stop')

  ##------------------------------------------ 

  def __should_continue(self,N_Iter_SelfCons_Max) :
    """ stop test"""
    should_continue = True
    if MPI.IS_MASTER_NODE():
      if (self.Iteration_Number > N_Iter_SelfCons_Max):
        should_continue = False
    should_continue = MPI.bcast(should_continue)
    return should_continue

  #**************************************************************************************************************

  
  def handler(self,signum, frame):
    print 'Signal handler called with signal', signum
    raise IOError, "Time elapsed"

  ##-------------------------------------------

  def run(self,N_Loops, Mixing_Coefficient = 0.5, MaxTime = 0 ):
    r"""
      Run the DMFT Loop with the following algorithm :: 
       
        while STOP_CONDITION : 
            self.Self_Consistency()
            for solver in Solver_List : S.solve()
            self.PostSolver() # defaults : does nothing

      where STOP_CONDITION is determined by the number of iterations.
        
      :param N_Loops:    Maximum number of iteration of the loop
      :param Mixing_Coefficient: 
      :param MaxTime: Maximum time of the loop.
    """

    # Set up the signal
    #   MPI.report("DMFTlab Job PID = %s"%os.getpid())
    # Set the signal handler and a 5-second alarm
    signal.signal(signal.SIGALRM, self.handler)
    signal.alarm(MaxTime)
 
    should_continue = True
    
    while (should_continue):
      MPI.report("------ Node : %d -------- Iteration Number = %d"%(MPI.rank,self.Iteration_Number))
      
      self.Self_Consistency()

      # call all solvers
      for n,sol in enumerate(self.SolverList) :
        if hasattr(self,"Chemical_potential") : sol.Chemical_potential=self.Chemical_potential
        sol.Iteration_Number=self.Iteration_Number
        sol.Solve()
        sol.Sigma  = sol.Sigma * Mixing_Coefficient + sol.Sigma_Old * (1-Mixing_Coefficient)
      
      # post-solver processing
      self.PostSolver()
                         
      self.Iteration_Number +=1
      should_continue = self.__should_continue(N_Loops)
 
    # end of the while loop
    MPI.report("----------- END of DMFT_Loop ----------------")
    MPI.barrier()
   

