from boost cimport *
from pytriqs.base.gf.local.gf cimport *
import pytriqs.base.gf.local
from types import *
from pytriqs.base.gf_local import *
from pytriqs.solvers import SolverBase
from pytriqs.solvers.operators import *
from pytriqs.base.utility.my_utils import *
import pytriqs.base.utility.parameters as parameters
import pytriqs.base.utility.mpi as mpi

cdef extern from "applications/impurity_solvers/ctqmc_hyb/Hloc.hpp":

    cdef cppclass Hloc_c "Hloc":
      Hloc_c(int, int, boost_object, boost_object, boost_object, boost_object, int)

cdef extern from "applications/impurity_solvers/ctqmc_hyb/MC.hpp" namespace "triqs::app::impurity_solvers":

    cdef cppclass solver_c "triqs::app::impurity_solvers::ctqmc_hyb":
      solver_c(boost_object, Hloc_c *, gf_block_imtime, gf_block_imtime, gf_block_imtime, gf_block_imtime, gf_block_legendre)
      void solve()


class Solver(SolverBase):
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
                "Random_Seed" : ("Seed for the random generator", 34788+928374*mpi.rank, IntType),
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

        parameters.check_no_parameters_not_in_union_of_dicts(param, self.Required, self.Optional)
        SolverBase.__init__(self,GFstruct,param)
        self.beta = float(Beta)
        self.Verbosity = 2 if mpi.rank ==0 else 0

        # Green function in frequencies
        a_list = [a for a,al in self.GFStruct]
        glist = [ GfImFreq(indices = al, beta = self.beta, n_matsubara =N_Matsubara_Frequencies) for a,al in self.GFStruct]
        self.G0 = BlockGf(name_list = a_list, block_list = glist, make_copies=False, name="G0")
        self.G = BlockGf(name_block_generator = self.G0, make_copies=True, name="G")
        self.F = BlockGf(name_block_generator = self.G0, make_copies=True, name="F")
        self.Sigma = BlockGf(name_block_generator = self.G0, make_copies=True, name="Sigma")
        self.Sigma_Old = BlockGf(name_block_generator = self.G0, make_copies=True, name="Sigma_Old")
        self.name = 'Hybridization Expansion'

         # first check that all indices of the Green Function do correspond to a C operator.
        for a,alpha_list in self.GFStruct :
          for alpha in alpha_list :
            if (a,alpha) not in operators.C_list_names() :
              raise "Error : Some indices (%s,%s) of the Green function do not correspond to existing operator"%(a,alpha)

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
        mpi.report(parameters.check(self.__dict__,self.Required,self.Optional))

        # We have to add the Hamiltonian the epsilon part of G0
        if type(self.H_Local) != type(Operator()) : raise "H_Local is not an operator"
        H = self.H_Local
        for a,alpha_list in  self.GFStruct :
            for mu in alpha_list : 
                for nu in alpha_list : 
                    H += real(self.G0[a]._tail[2][mu,nu]) * Cdag(a,mu)*C(a,nu)

        OPdict = {"Hamiltonian": H}
        mpi.report("Hamiltonian with Eps0 term  : ",H)
        
        # First separate the quantum Numbers that are operators and those which are symmetries.
        QuantumNumberOperators  = dict( (n,op) for (n,op) in self.Quantum_Numbers.items() if type(op) == type(Operator()))
        QuantumNumberSymmetries = dict( (n,op) for (n,op) in self.Quantum_Numbers.items() if type(op) != type(Operator()))

        # Check that the quantum numbers commutes with the Hamiltonian
        for name,op in QuantumNumberOperators.items():
            assert commutator(self.H_Local ,op).is_zero(), "One quantum number is not commuting with Hamiltonian"
            OPdict[name]=op

        # Complete the OPdict with the fundamental operators
        OPdict, nf, nb, SymChar, NameOpFundamentalList = operators.complete_op_list_with_fundamentals(OPdict)

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
              mpi.report("Operator %s already defined as %s, using this instead for measuring"%(name,opn))
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
              mpi.report("Operator %s already defined as %s, using this instead for measuring"%(name,opn))
              if opn not in self.OpCorr_To_Average_List: self.OpCorr_To_Average_List.append(opn)
        # Create storage for data:
        Nops = len(self.OpCorr_To_Average_List)
        f = lambda L : GfImTime(indices = [0], beta = self.beta, n_time_points =L )
        if (Nops>0):
            self.Measured_Time_Correlators_Results = BlockGf(name_block_generator = [ ( n,f(self.Measured_Time_Correlators[n][1]) ) for n in self.Measured_Time_Correlators], make_copies=False)
        else:
            self.Measured_Time_Correlators_Results = BlockGf(name_block_generator = [ ( 'OpCorr',f(2) ) ], make_copies=False)

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
                F = extend_function_on_fundamentals(Map_GM_to_Fund_Ops(GM))
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
            F = extend_function_on_fundamentals(Map_GM_to_Fund_Ops(GM))
            m = {}
            for name,op in OPdict.items() :
                op_im = F(op)
                n1,n2 = myfind(op),myfind(op_im)
                m[n1] = n2
            name = "%s"%n
            self.Global_Moves_Mapping_List.append((proba,m,name))
        #mpi.report ("Global_Moves_Mapping_List", self.Global_Moves_Mapping_List)

        # Now add the operator for F calculation if needed
        if self.Use_F :
            Hloc_WithoutQuadratic = self.H_Local.remove_quadratic()
            for n,op in OPdict.items() :
                if op.is_Fundamental():
                    op2 = commutator(Hloc_WithoutQuadratic,op)
                    if not mysearch(op2) : OPdict["%s_Comm_Hloc"%n] = op2

        # All operators have real coefficients. Check this and remove the 0j term
        # since the C++ expects operators with real numbers 
        for n,op in OPdict.items(): op.make_coef_real_and_check()

        # Transcription of operators for C++
        Oplist2 = operators.transcribe_op_list_for_C(OPdict)
        SymList = [sym for (n,sym) in SymChar.items() if n in QuantumNumberSymmetries]
        #self.H_diag = C_Module.Hloc(nf,nb,Oplist2,QuantumNumberOperators,SymList,self.Quantum_Numbers_Selection,0) 

        # Create the C_Cag_Ops array which describes the grouping of (C,Cdagger) operator
        # for the MonteCarlo moves : (a, alpha) block structure [ [ (C_name, Cdag_name)]]
        self.C_Cdag_Ops = [ [ (myfind(C(a,alpha)), myfind(Cdag(a,alpha))) for alpha in al ] for a,al in self.GFStruct]

        # Define G0_inv and correct it to have G0 to have perfect 1/omega behavior
        self.G0_inv = inverse(self.G0)
        Delta = self.G0_inv.delta()
        for n,g in self.G0_inv:
          assert(g.N1==g.N2)
          identity=numpy.identity(g.N1)
          self.G0[n] <<= gf_init.A_Omega_Plus_B(identity, g._tail[0])
          self.G0[n] -= Delta[n]
          #self.G0[n] <<= iOmega_n + g._tail[0] - Delta[n]
        self.G0_inv <<= self.G0
        self.G0.invert()

        # Construct the function in tau
        f = lambda g,L : GfImTime(indices = g.indices, beta = g.beta, n_time_points =L )
        self.Delta_tau = BlockGf(name_block_generator = [ (n,f(g,self.N_Time_Slices_Delta) )   for n,g in self.G], make_copies=False, name='D')
        self.G_tau = BlockGf(name_block_generator = [ (n,f(g,self.N_Time_Slices_Gtau) )    for n,g in self.G], make_copies=False, name='G')
        self.F_tau = BlockGf(name_block_generator = self.G_tau, make_copies=True, name='F')
        
        for (i,gt) in self.Delta_tau : gt.set_from_inverse_fourier(Delta[i])
        mpi.report("Inv Fourier done")
        if (self.Legendre_Accumulation):
            self.G_Legendre = BlockGf(name_block_generator = [ (n,GfLegendre(indices =g.indices, beta =g.beta, n_legendre_coeffs =self.N_Legendre_Coeffs) )   for n,g in self.G], make_copies=False, name='Gl')
        else:
            self.G_Legendre = BlockGf(name_block_generator = [ (n,GfLegendre(indices =[1], beta =g.beta, n_legendre_coeffs =1) ) for n,g in self.G], make_copies=False, name='Gl') # G_Legendre must not be empty but is not needed in this case. So I make it as small as possible.
        
        # Starting the C++ code
        self.Sigma_Old <<= self.Sigma

        # C++ solver
        f = lambda g, L: pytriqs.base.gf.local.GfImTime(indices=g.indices, beta=g.beta, n_time_points=L)
        DD = pytriqs.base.gf.local.BlockGf(name_block_generator = [ (n,f(g,self.N_Time_Slices_Delta)) for n,g in self.Delta_tau], make_copies=False)
        GG = pytriqs.base.gf.local.BlockGf(name_block_generator = [ (n,f(g,self.N_Time_Slices_Gtau)) for n,g in self.G_tau], make_copies=False)
        FF = pytriqs.base.gf.local.BlockGf(name_block_generator = GG, make_copies=True, name='F')

        f = lambda L : pytriqs.base.gf.local.GfImTime(indices = [0], beta = self.beta, n_time_points =L )
        if (Nops>0):
          MM = pytriqs.base.gf.local.BlockGf(name_block_generator = [ (n,f(self.Measured_Time_Correlators[n][1]) ) for n in self.Measured_Time_Correlators], make_copies=False)
        else:
          MM = pytriqs.base.gf.local.BlockGf(name_block_generator = [ ('OpCorr',f(2)) ], make_copies=False)

        f = lambda g, L: pytriqs.base.gf.local.GfLegendre(indices=g.indices, beta=g.beta, n_legendre_points=L)
        LL = pytriqs.base.gf.local.BlockGf(name_block_generator = [ (n,f(g,self.N_Legendre_Coeffs)) for n,g in self.G_Legendre], make_copies=False)

        for n,d in DD: d.data[:,:,:] = self.Delta_tau[n]._data.array[:,:,:]

        solver_c(boost_object(self.__dict__),
                 new Hloc_c(nf, nb, boost_object(Oplist2), boost_object(QuantumNumberOperators), boost_object(SymList),
                            boost_object(self.Quantum_Numbers_Selection), 0),
                 as_gf_block_imtime(GG),
                 as_gf_block_imtime(FF),
                 as_gf_block_imtime(DD),
                 as_gf_block_imtime(MM),
                 as_gf_block_legendre(LL)).solve()

        for n,g in self.G_Legendre: g._data.array[:,:,:] = LL[n].data[:,:,:]
        for n,g in self.G_tau: g._data.array[:,:,:] = GG[n].data[:,:,:]
        for n,g in self.F_tau: g._data.array[:,:,:] = FF[n].data[:,:,:]
        for n,g in self.Measured_Time_Correlators_Results: g._data.array[:,:,:] = MM[n].data[:,:,:]
        
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
              self.G[name].set_from_fourier(g)

          # This is very sick... but what can we do???
          self.Sigma <<= self.G0_inv - inverse(self.G)
          self.fitTails()
          self.G <<= inverse(self.G0_inv - self.Sigma)

        # Now find the self-energy
        self.Sigma <<= self.G0_inv - inverse(self.G)

        mpi.report("Solver %(name)s has ended."%self.__dict__)

        # for operator averages: if twice defined operator, rename output:
        for op1,op2 in self.twice_defined_Ops.items():
            self.Measured_Operators_Results[op1] = self.Measured_Operators_Results[op2]
        for op1,op2 in self.twice_defined_Ops.items():
            if op2 in self.Measured_Operators_Results.keys(): del self.Measured_Operators_Results[op2]

        if self.Use_F :
            for (n,f) in self.F: f.set_from_fourier(self.F_tau[n])
            self.G2 = self.G0 + self.G0 * self.F
            self.Sigma2 = self.F * inverse(self.G2)
            
    #--------------------------------------------------

    def fitTails(self):
        # fits the tails of the noise monte carlo data

        for n,sig in self.Sigma:

            known_coeff = numpy.zeros([sig.N1,sig.N2,1],numpy.float_)
            msh = [x.imag for x in sig.mesh]
            fit_start = msh[self.fitting_Frequency_Start]
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
      if issue_warning: mpi.report(msg)
