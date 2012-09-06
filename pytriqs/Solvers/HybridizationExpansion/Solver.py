
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

from types import *
from pytriqs.Base.GF_Local import *
from pytriqs.Solvers import Solver_Base
from pytriqs.Solvers.Operators import *
from pytriqs.Base.Utility.myUtils import *
import pytriqs.Base.Utility.Parameters as Parameters
import pytriqs.Base.Utility.MPI as MPI

try :
    import pytriqs_Solver_HybridizationExpansion as C_Module
except : 
    raise ImportError, "Oops ! Did you compile the Hybridization solver ? I can't find it !"

__add__ = []

class Solver(Solver_Base):
    """
    Hybridization QMC solver.
    """

    Required = {
                "N_Cycles" : ("Maximum number of Monte Carlo cycles", IntType),
                "H_Local" : ("Local Hamiltonian", InstanceType),
                "Quantum_Numbers" : ("Quantum Numbers", DictType),
                }

    Optional = {
                "Random_Generator_Name" : ("Name of the random number generator", "", StringType),
                "Random_Seed" : ("Seed for the random generator", 34788+928374*MPI.rank, IntType),
                "Length_Cycle":("Length of one QMC cycle", 200, IntType),
                "Legendre_Accumulation" : ("Do we accumulate in legendre?", True, BooleanType),
                "N_Legendre_Coeffs" : ("Number of Legendre coefficients that are used in practice", 50, IntType),
                "Time_Accumulation" : ("Do we accumulate in imaginary-time?", False, BooleanType),
                "N_Warmup_Cycles" : ("Number of warming iterations (cycles)", 1000, IntType),
                "N_Frequencies_Accumulated" : ("Number of frequencies to accumulate", 100, IntType),
                "Fitting_Frequency_Start" : ("Frequency at which the fit starts", 50, IntType),
                "N_Time_Slices_Delta" : ("Number of times slices in Delta", 10000, IntType),
                "N_Time_Slices_Gtau" : ("Number of times slices in G_tau", 10000, IntType),
                "Nmax_Matrix" : ("Initial size of the determinant matrices", 100, IntType),
                "Global_Moves" : ("A list of global moves", [], ListType),
                "Use_Segment_Picture" : ("Guarantee is made that the C^+ C alternate in each config", False, BooleanType),
                "Proba_Insert_Remove" : ("Probability to insert/remove operators", 1.0, FloatType),
                "Proba_Move" : ("Probability to move operators", 1.0, FloatType),
                "Measured_Operators" : ("A dict of operators that will be averaged", {}, DictType),
                "Measured_Time_Correlators" : ("A dict of operators, whose time correlations are to be measured", {}, DictType),
                "Record_Statistics_Configurations" : ("(Expert only) Get the kink length statistics", False, BooleanType),
                "Keep_Full_MC_Series" : ("(Expert only) Store the Green's function for later analysis", False, BooleanType),
                "Use_F" : ("(Expert only) Compute F", False, BooleanType),
                "Eta": ("(Expert only) Value of eta, the minimum value of the det", 0.0, FloatType),
                "Quantum_Numbers_Selection" : ("(Prototype) A function to select quantum numbers", lambda qn: True, FunctionType),
                }

    Deprecated = ()

    Checks =  ()

    #--------------------------------------------------

    def __init__(self,Beta,GFstruct,N_Matsubara_Frequencies=1025,**param):
        """
        :param Beta: The inverse temperature
        :param GFstruct: The structure of the Green's functions. It must be a list of tuples,
                         each representing a block of the Green's function. The tuples have two
                         elements, the name of the block and a list of indices. For example:
                         [ ('up', [1,2,3]),  ('down', [1,2,3]) ].
        :param N_Matsubara_Frequencies: (Optional, default = 1025) How many Matsubara frequencies
                                        are used for the Green's functions.
        :param param: A list of extra parameters, described below.
        """

        # For backward compatibility
        self.update_params(param)

        Parameters.check_no_parameters_not_in_union_of_dicts(param, self.Required, self.Optional)
        Solver_Base.__init__(self,GFstruct,param)
        self.Beta = float(Beta)
        self.Verbosity = 2 if MPI.rank ==0 else 0

        # Green function in frequencies
        a_list = [a for a,al in self.GFStruct]
        glist = [ GFBloc_ImFreq(Indices = al, Beta = self.Beta, NFreqMatsubara=N_Matsubara_Frequencies) for a,al in self.GFStruct]
        self.G0 = GF(NameList = a_list, BlockList = glist, Copy=False, Name="G0")
        self.G = GF(Name_Block_Generator = self.G0, Copy=True, Name="G")
        self.F = GF(Name_Block_Generator = self.G0, Copy=True, Name="F")
        self.Sigma = GF(Name_Block_Generator = self.G0, Copy=True, Name="Sigma")
        self.Sigma_Old = GF(Name_Block_Generator = self.G0, Copy=True, Name="Sigma_Old")
        self.Name = 'Hybridization Expansion'

         # first check that all indices of the Green Function do correspond to a C operator.
        for a,alpha_list in self.GFStruct :
          for alpha in alpha_list :
            if (a,alpha) not in Operators.ClistNames() :
              raise "Error : Some indices (%s,%s) of the Green function do not correspond to existing operator"%(a,alpha)

    #--------------------------------------------------

    @classmethod
    def Random_Generators_Available(cls) :
        """ List of all random number generators compiled """
        return C_Module.Random_Generators_Available(' ').split()

    #--------------------------------------------------
                
    def Solve(self):
        """ Solve the impurity problem """

        # Find if an operator is in oplist
        def mysearch(op):
            l = [ k for (k,v) in OPdict.items() if (v-op).is_zero()]
            assert len(l) <=1
            return l[0] if l else None

        # Same but raises an error if pb
        def myfind(op):
            r = mysearch(op)
            if r==None : raise "Operator %s can not be found by myfind !"%r
            return r

        # For backward compatibility
        self.update_params(self.__dict__)

        # Test all a parameters before solutions
        MPI.report(Parameters.check(self.__dict__,self.Required,self.Optional))

        # We have to add the Hamiltonian the epsilon part of G0
        if type(self.H_Local) != type(Operator()) : raise "H_Local is not an operator"
        H = self.H_Local
        for a,alpha_list in  self.GFStruct :
            for mu in alpha_list : 
                for nu in alpha_list : 
                    H += real(self.G0[a]._tail[2][mu,nu]) * Cdag(a,mu)*C(a,nu)

        OPdict = {"Hamiltonian": H}
        MPI.report("Hamiltonian with Eps0 term  : ",H)
        
        # First separate the quantum Numbers that are operators and those which are symmetries.
        QuantumNumberOperators  = dict( (n,op) for (n,op) in self.Quantum_Numbers.items() if type(op) == type(Operator()))
        QuantumNumberSymmetries = dict( (n,op) for (n,op) in self.Quantum_Numbers.items() if type(op) != type(Operator()))

        # Check that the quantum numbers commutes with the Hamiltonian
        for name,op in QuantumNumberOperators.items():
            assert Commutator(self.H_Local ,op).is_zero(), "One quantum number is not commuting with Hamiltonian"
            OPdict[name]=op

        # Complete the OPdict with the fundamental operators
        OPdict, nf, nb, SymChar, NameOpFundamentalList = Operators.Complete_OperatorsList_with_Fundamentals(OPdict)

        # Add the operators to be averaged in OPdict and prepare the list for the C-code
        self.Measured_Operators_Results = {}
        self.twice_defined_Ops = {}
        self.Operators_To_Average_List = []
        for name, op in self.Measured_Operators.items():
          opn = mysearch(op)
          if opn == None : 
              OPdict[name] = op
              self.Measured_Operators_Results[name] = 0.0
              self.Operators_To_Average_List.append(name)
          else:
              MPI.report("Operator %s already defined as %s, using this instead for measuring"%(name,opn))
              self.twice_defined_Ops[name] = opn
              self.Measured_Operators_Results[opn] = 0.0
              if opn not in self.Operators_To_Average_List: self.Operators_To_Average_List.append(opn)

        # Time correlation functions are added
        self.OpCorr_To_Average_List = []
        for name, op in self.Measured_Time_Correlators.items():
          opn = mysearch(op[0])
          if opn == None : 
              OPdict[name] = op[0]
              self.OpCorr_To_Average_List.append(name)
          else:
              MPI.report("Operator %s already defined as %s, using this instead for measuring"%(name,opn))
              if opn not in self.OpCorr_To_Average_List: self.OpCorr_To_Average_List.append(opn)
        # Create storage for data:
        Nops = len(self.OpCorr_To_Average_List)
        f = lambda L : GFBloc_ImTime(Indices= [0], Beta = self.Beta, NTimeSlices=L )
        if (Nops>0):
            self.Measured_Time_Correlators_Results = GF(Name_Block_Generator = [ ( n,f(self.Measured_Time_Correlators[n][1]) ) for n in self.Measured_Time_Correlators], Copy=False)
        else:
            self.Measured_Time_Correlators_Results = GF(Name_Block_Generator = [ ( 'OpCorr',f(2) ) ], Copy=False)

        # Take care of the global moves

        # First, given a function (a,alpha,dagger) -> (a', alpha', dagger')
        # I construct a function on fundamental operators
        def Map_GM_to_Fund_Ops( GM ) :
            def f(fop) :
                a,alpha, dagger = fop.name + (fop.dag,)
                ap,alphap,daggerp = GM((a,alpha,dagger))
                return Cdag(ap,alphap) if daggerp else C(ap,alphap)
            return f

        # Complete the OpList so that it is closed under the global moves
        while 1:
            added_something = False
            for n,(proba,GM) in enumerate(self.Global_Moves):
                # F is a function that map all operators according to the global move
                F = Extend_Function_on_Fundamentals(Map_GM_to_Fund_Ops(GM))
                # Make sure that OPdict is complete, i.e. all images of OPdict operators are in OPdict
                for name,op in OPdict.items() :
                    op_im = F(op)
                    if mysearch(op_im)==None :
                        # find the key and put in in the dictionnary
                        i=0
                        while 1:
                            new_name = name + 'GM' +  i*'_' + "%s"%n
                            if new_name not in OPdict : break
                        added_something = True
                        OPdict[new_name] = op_im
            # break the while loop
            if not added_something: break

        # Now I have all operators, I make the transcription of the global moves
        self.Global_Moves_Mapping_List = []
        for n,(proba,GM) in enumerate(self.Global_Moves):
            F = Extend_Function_on_Fundamentals(Map_GM_to_Fund_Ops(GM))
            m = {}
            for name,op in OPdict.items() :
                op_im = F(op)
                n1,n2 = myfind(op),myfind(op_im)
                m[n1] = n2
            name = "%s"%n
            self.Global_Moves_Mapping_List.append((proba,m,name))
        #MPI.report ("Global_Moves_Mapping_List", self.Global_Moves_Mapping_List)

        # Now add the operator for F calculation if needed
        if self.Use_F :
            Hloc_WithoutQuadratic = self.H_Local.RemoveQuadraticTerms()
            for n,op in OPdict.items() :
                if op.is_Fundamental():
                    op2 = Commutator(Hloc_WithoutQuadratic,op)
                    if not mysearch(op2) : OPdict["%s_Comm_Hloc"%n] = op2

        # All operators have real coefficients. Check this and remove the 0j term
        # since the C++ expects operators with real numbers 
        for n,op in OPdict.items(): op.make_coef_real_and_check()

        # Transcription of operators for C++
        Oplist2 = Operators.Transcribe_OpList_for_C(OPdict)
        SymList = [sym for (n,sym) in SymChar.items() if n in QuantumNumberSymmetries]
        self.H_diag = C_Module.Hloc(nf,nb,Oplist2,QuantumNumberOperators,SymList,self.Quantum_Numbers_Selection,0) 

        # Create the C_Cag_Ops array which describes the grouping of (C,Cdagger) operator
        # for the MonteCarlo moves : (a, alpha) block structure [ [ (C_name, Cdag_name)]]
        self.C_Cdag_Ops = [ [ (myfind(C(a,alpha)), myfind(Cdag(a,alpha))) for alpha in al ] for a,al in self.GFStruct]

        # Define G0_inv and correct it to have G0 to have perfect 1/omega behavior
        self.G0_inv = inverse(self.G0)
        Delta = self.G0_inv.Delta()
        for n,g in self.G0_inv:
          assert(g.N1==g.N2)
          identity=numpy.identity(g.N1)
          self.G0[n] <<= GF_Initializers.A_Omega_Plus_B(identity, g._tail[0])
          self.G0[n] -= Delta[n]
          #self.G0[n] <<= iOmega_n + g._tail[0] - Delta[n]
        self.G0_inv <<= self.G0
        self.G0.invert()

        # Construct the function in tau
        f = lambda g,L : GFBloc_ImTime(Indices= g.Indices, Beta = g.Beta, NTimeSlices=L )
        self.Delta_tau = GF(Name_Block_Generator = [ (n,f(g,self.N_Time_Slices_Delta) )   for n,g in self.G], Copy=False, Name='D')
        self.G_tau = GF(Name_Block_Generator = [ (n,f(g,self.N_Time_Slices_Gtau) )    for n,g in self.G], Copy=False, Name='G')
        self.F_tau = GF(Name_Block_Generator = self.G_tau, Copy=True, Name='F')
        
        for (i,gt) in self.Delta_tau : gt.setFromInverseFourierOf(Delta[i])
        MPI.report("Inv Fourier done")
        if (self.Legendre_Accumulation):
            self.G_Legendre = GF(Name_Block_Generator = [ (n,GFBloc_ImLegendre(Indices=g.Indices, Beta=g.Beta, NLegendreCoeffs=self.N_Legendre_Coeffs) )   for n,g in self.G], Copy=False, Name='Gl')
        else:
            self.G_Legendre = GF(Name_Block_Generator = [ (n,GFBloc_ImLegendre(Indices=[1], Beta=g.Beta, NLegendreCoeffs=1) ) for n,g in self.G], Copy=False, Name='Gl') # G_Legendre must not be empty but is not needed in this case. So I make it as small as possible.
        
        # Starting the C++ code
        self.Sigma_Old <<= self.Sigma
        C_Module.MC_solve(self.__dict__ ) # C++ solver
        
        # Compute G on Matsubara axis possibly fitting the tail
        if self.Legendre_Accumulation:
          for s,g in self.G:
            identity=numpy.zeros([g.N1,g.N2],numpy.float)
            for i,m in enumerate (g._IndicesL):
              for j,n in enumerate (g._IndicesR):
                if m==n: identity[i,j]=1
            self.G_Legendre[s].enforce_discontinuity(identity) # set the known tail
            g <<= LegendreToMatsubara(self.G_Legendre[s])
        else:
          if (self.Time_Accumulation):
            for name, g in self.G_tau:
              identity=numpy.zeros([g.N1,g.N2],numpy.float)
              for i,m in enumerate (g._IndicesL):
                for j,n in enumerate (g._IndicesR):
                  if m==n: identity[i,j]=1
              g._tail.zero()
              g._tail[1] = identity
              self.G[name].setFromFourierOf(g)

          # This is very sick... but what can we do???
          self.Sigma <<= self.G0_inv - inverse(self.G)
          self.fitTails()
          self.G <<= inverse(self.G0_inv - self.Sigma)

        # Now find the self-energy
        self.Sigma <<= self.G0_inv - inverse(self.G)

        MPI.report("Solver %(Name)s has ended."%self.__dict__)

        # for operator averages: if twice defined operator, rename output:
        for op1,op2 in self.twice_defined_Ops.items():
            self.Measured_Operators_Results[op1] = self.Measured_Operators_Results[op2]
        for op1,op2 in self.twice_defined_Ops.items():
            if op2 in self.Measured_Operators_Results.keys(): del self.Measured_Operators_Results[op2]

        if self.Use_F :
            for (n,f) in self.F: f.setFromFourierOf(self.F_tau[n])
            self.G2 = self.G0 + self.G0 * self.F
            self.Sigma2 = self.F * inverse(self.G2)
            
    #--------------------------------------------------

    def fitTails(self):
        # fits the tails of the noise monte carlo data

        for n,sig in self.Sigma:

            known_coeff = numpy.zeros([sig.N1,sig.N2,1],numpy.float_)
            msh = [x.imag for x in sig.mesh]
            fit_start = msh[self.Fitting_Frequency_Start]
            fit_stop  = msh[self.N_Frequencies_Accumulated-1]
            
            sig.fitTail(fixed_coef = known_coeff, order_max = 3, fit_start = fit_start, fit_stop = fit_stop)

    #--------------------------------------------------

    # This is for backward compatibility
    def update_params(self, d):

      allparams = [
        ('QMC_N_cycles_MAX', 'N_Cycles'),
        ('NCycles', 'N_Cycles'),
        ('Hloc', 'H_Local'),
        ('QuantumNumbers', 'Quantum_Numbers'),
        ('Length_One_QMC_Cycle', 'Length_Cycle'),
        ('Number_Warming_Iteration', 'N_Warmup_Cycles'),
        ('Number_Frequencies_Accumulated', 'N_Frequencies_Accumulated'),
        ('Global_Move', 'Global_Moves'),
        ('UseSegmentPicture', 'Use_Segment_Picture'),
        ('Proba_Move_Insert_Remove_Kink', 'Proba_Insert_Remove'),
        ('Proba_Move_Move_Kink', 'Proba_Move'),
        ('OperatorsToAverage', 'Measured_Operators'),
        ('OpCorrToAverage', 'Measured_Time_Correlators'),
        ('KeepGF_MC_series', 'Keep_Full_MC_Series'),
        ('DecorrelationAnalysisG_NFreq', 'Decorrelation_Analysis_G_NFreq'),
        ('RecordStatisticConfigurations', 'Record_Statistics_Configurations')
      ]

      issue_warning = False
      for (old, new) in allparams:
        if old in d:
          val = d.pop(old)
          d.update({new:val})
          issue_warning = True

      msg = """
**********************************************************************************
 Warning: some parameters you used have a different name now and will be
 deprecated in future versions. Please check the documentation.
**********************************************************************************
"""
      if issue_warning: MPI.report(msg)
