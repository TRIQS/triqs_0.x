
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

__all__ = ["DMFTLoopGeneric"]

import shelve,os,signal,types
from pytriqs.applications.impurity_solvers import *
import pytriqs.utility.mpi as mpi

class DMFTLoopGeneric:
  """ This class provides a generic and simple DMFT loop for convenience.
      Given a list of solvers, its method `run` runs the generic DMFT loop.

      Usage : 

        * Derive and redefine Self_Consistency() (and maybe post_solver)

      Useful elements : 

        * self.Iteration_Number is the number of the DMFT loop iteration
              
      """

  def __init__(self, solver_list,  **params) : 
    """

    :param solver_list:  * List of the solvers to be called in the loop. 
                         *  If there is only one solver S, the syntax solver_list = S 
                            is tolerated instead of solver_list = [S]

    :param params: 
                        * all keywords arguments are simply put in the namespace
                        * of the loop e.g. :: 
                            
                             myloop(solver_list = S,Chemical_potential= 2.0, theta = 78)
                          
                          myloop will be constructed with ::

                            self.Chemical_potential = 2.0
                            self.theta = 78

    """
    import warnings
    warnings.warn("The class %s is deprecated and scheduled for removal in later version of TRIQS. Please upgrade your code (write your own loop)"%(self.__class__.__name__), DeprecationWarning)
    self.SolverList = solver_list if type(solver_list) in [types.ListType,types.TupleType] else [solver_list]
    for s in self.SolverList : assert isinstance(s,SolverBase), "Oops, one element of solver_list is not a solver !!!"

    self.Iteration_Number  = 1
    # Put the additionnal parameters in my name space.
    self.__dict__.update(params)
     
  ##--------------------------------------------
      
  def post_solver(self) : 
      """ Do nothing by default """
      pass
     
  ##-------------------------------------------

  # Obsolete : due to old RUPC config.
  def soft_kill_signal(self) : return os.path.exists('.stop')

  ##------------------------------------------ 

  def __should_continue(self, n_iter_max) :
    """ stop test"""
    should_continue = True
    if mpi.is_master_node():
      if (self.Iteration_Number > n_iter_max):
        should_continue = False
    should_continue = mpi.bcast(should_continue)
    return should_continue

  #**************************************************************************************************************

  
  def handler(self,signum, frame):
    print 'Signal handler called with signal', signum
    raise IOError, "Time elapsed"

  ##-------------------------------------------

  def run(self, n_loops, mixing_coeff = 0.5, max_time = 0 ):
    r"""
      Run the DMFT Loop with the following algorithm :: 
       
        while STOP_CONDITION : 
            self.Self_Consistency()
            for solver in solver_list : S.solve()
            self.post_solver() # defaults : does nothing

      where STOP_CONDITION is determined by the number of iterations.
        
      :param n_loops:    Maximum number of iteration of the loop
      :param mixing_coeff: 
      :param max_time: Maximum time of the loop.
    """

    # Set up the signal
    #   mpi.report("DMFTlab Job PID = %s"%os.getpid())
    # Set the signal handler and a 5-second alarm
    signal.signal(signal.SIGALRM, self.handler)
    signal.alarm(max_time)
 
    should_continue = True
    
    while (should_continue):
      mpi.report("------ Node : %d -------- Iteration Number = %d"%(mpi.rank,self.Iteration_Number))
      
      self.Self_Consistency()

      # call all solvers
      for n,sol in enumerate(self.SolverList) :
        if hasattr(self,"Chemical_potential") : sol.Chemical_potential=self.Chemical_potential
        sol.Iteration_Number=self.Iteration_Number
        sol.Solve()
        sol.Sigma  = sol.Sigma * mixing_coeff + sol.Sigma_Old * (1-mixing_coeff)
      
      # post-solver processing
      self.post_solver()
                         
      self.Iteration_Number +=1
      should_continue = self.__should_continue(n_loops)
 
    # end of the while loop
    mpi.report("----------- END of DMFT_Loop ----------------")
    mpi.barrier()
   

