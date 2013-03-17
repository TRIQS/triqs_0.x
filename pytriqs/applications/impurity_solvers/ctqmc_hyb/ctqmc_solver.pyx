from boost cimport *
from pytriqs.gf.local.gf cimport *
from types import *
from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import A_Omega_Plus_B
from pytriqs.applications.impurity_solvers.operators import *
import pytriqs.utility.mpi as mpi
import numpy

cdef extern from "applications/impurity_solvers/ctqmc_hyb/hloc.hpp":

    cdef cppclass Hloc_c "Hloc":
      Hloc_c(int, int, boost_object, boost_object, boost_object, boost_object, int)

cdef extern from "applications/impurity_solvers/ctqmc_hyb/ctqmc.hpp" namespace "triqs::app::impurity_solvers":

    cdef cppclass solver_c "triqs::app::impurity_solvers::ctqmc_hyb":
      solver_c(boost_object, Hloc_c *, gf_block_imtime, gf_block_imtime, gf_block_imtime, gf_block_imtime, gf_block_legendre)
      void solve()


class Solver:
    """
    Hybridization-expansion QMC solver.
    """

    def __init__(self, beta, gf_struct, n_matsubara=1025):
        """
        :param beta: The inverse temperature
        :param gf_struct: The structure of the Green's functions. It must be a list of tuples,
                         each representing a block of the Green's function. The tuples have two
                         elements, the name of the block and a list of indices. For example:
                         [ ('up', [1,2,3]),  ('down', [1,2,3]) ].
        :param n_matsubara: (Optional, default = 1025) How many Matsubara frequencies
                                        are used for the Green's functions.
        """

        # Set some useful variables
        self.beta = float(beta)
        self.gf_struct = gf_struct[:]
        self.name = 'Hybridization Expansion'

        # Green's functions in frequencies
        a_list = [a for a,al in self.gf_struct]
        glist = [ GfImFreq(indices = al, beta = self.beta, n_points = n_matsubara) for a,al in self.gf_struct]
        self.G0 = BlockGf(name_list = a_list, block_list = glist, make_copies=False, name="G0")
        self.G = BlockGf(name_block_generator = self.G0, make_copies=True, name="G")
        self.F = BlockGf(name_block_generator = self.G0, make_copies=True, name="F")
        self.Sigma = BlockGf(name_block_generator = self.G0, make_copies=True, name="Sigma")
        self.Sigma_old = BlockGf(name_block_generator = self.G0, make_copies=True, name="Sigma_old")

    #--------------------------------------------------

    def solve(self, **params):
        """ Solve the impurity problem """

        # Required parameters
        self.H_local = params.pop("H_local")
        self.quantum_numbers = params.pop("quantum_numbers")

        # Monte Carlo part
        self.n_cycles = params.pop("n_cycles")
        self.length_cycle = params.pop("length_cycle", 200)
        self.n_warmup_cycles = params.pop("n_warmup_cycles", 1000)
        self.random_name = params.pop("random_name", "")
        self.random_seed = params.pop("random_seed", 34788+928374*mpi.rank)
        self.verbosity = 2 if mpi.rank ==0 else 0
        self.max_time = params.pop("max_time",-1)

        # Other control parameters
        self.use_segment_picture = params.pop("use_segment_picture", False)
        self.legendre_accumulation = params.pop("legendre_accumulation", True)
        self.n_legendre = params.pop("n_legendre", 50)
        self.time_accumulation = params.pop("time_accumulation", False)
        self.fit_start = params.pop("fit_start", 50)
        self.fit_stop = params.pop("fit_stop", 100)
        self.n_time_slices_delta = params.pop("n_time_slices_delta", 10000)
        self.n_time_slices_gtau = params.pop("n_time_slices_gtau", 10000)
        self.n_max_matrix = params.pop("n_max_matrix", 100)
        self.global_moves = params.pop("global_moves", [])
        self.prob_insert_remove = params.pop("prob_insert_remove", 1.0)
        self.prob_move = params.pop("prob_move", 1.0)
        self.measured_operators = params.pop("measured_operators", {})
        self.measured_time_correlators = params.pop("measured_time_correlators", {})
        self.record_stat = params.pop("record_stat", False)
        self.use_f = params.pop("use_f", False)
        self.quantum_number_selection = params.pop("quantum_number_selection", lambda qn: True)

        # Check params has been consumed
        if params:
          print "\nWARNING: Some parameters have not been used:"
          for p in params: print "--> %s"%p
          print "\n"

        # Check that all indices of the Green's Function do correspond to a C operator
        for a,alpha_list in self.gf_struct :
          for alpha in alpha_list :
            if (a,alpha) not in operators.C_list_names() :
              print (a,alpha)
              raise "Error : Some indices (%s,%s) of the Green function do not correspond to existing operator"%(a,alpha)

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

        # We have to add the Hamiltonian the epsilon part of G0
        if type(self.H_local) != type(Operator()) : raise "H_local is not an operator"
        H = self.H_local
        for a,alpha_list in  self.gf_struct :
            for mu, amu in enumerate(alpha_list) :
                for nu, anu in enumerate(alpha_list) :
                    H += (self.G0[a].tail[2][mu,nu]).real * Cdag(a,amu)*C(a,anu)

        OPdict = {"Hamiltonian": H}
        mpi.report("Hamiltonian with Eps0 term  : ",H)

        # First separate the quantum Numbers that are operators and those which are symmetries.
        QuantumNumberOperators  = dict( (n,op) for (n,op) in self.quantum_numbers.items() if type(op) == type(Operator()))
        QuantumNumberSymmetries = dict( (n,op) for (n,op) in self.quantum_numbers.items() if type(op) != type(Operator()))

        # Check that the quantum numbers commutes with the Hamiltonian
        for name,op in QuantumNumberOperators.items():
            assert commutator(self.H_local ,op).is_zero(), "One quantum number is not commuting with Hamiltonian"
            OPdict[name]=op

        # Complete the OPdict with the fundamental operators
        OPdict, nf, nb, SymChar, NameOpFundamentalList = operators.complete_op_list_with_fundamentals(OPdict)

        # Add the operators to be averaged in OPdict and prepare the list for the C-code
        self.measured_operators_results = {}
        self.twice_defined_ops = {}
        self.operators_av_list = []
        for name, op in self.measured_operators.items():
          opn = mysearch(op)
          if opn == None :
              OPdict[name] = op
              self.measured_operators_results[name] = 0.0
              self.operators_av_list.append(name)
          else:
              mpi.report("Operator %s already defined as %s, using this instead for measuring"%(name,opn))
              self.twice_defined_ops[name] = opn
              self.measured_operators_results[opn] = 0.0
              if opn not in self.operators_av_list: self.operators_av_list.append(opn)

        # Time correlation functions are added
        self.opcorr_av_list = []
        for name, op in self.measured_time_correlators.items():
          opn = mysearch(op[0])
          if opn == None :
              OPdict[name] = op[0]
              self.opcorr_av_list.append(name)
          else:
              mpi.report("Operator %s already defined as %s, using this instead for measuring"%(name,opn))
              if opn not in self.opcorr_av_list: self.opcorr_av_list.append(opn)
        # Create storage for data:
        Nops = len(self.opcorr_av_list)
        f = lambda L : GfImTime(indices = [0], beta = self.beta, n_points =L )
        if (Nops>0):
            self.measured_time_correlators_results = BlockGf(name_block_generator = [ ( n,f(self.measured_time_correlators[n][1]) ) for n in self.measured_time_correlators], make_copies=False)
        else:
            self.measured_time_correlators_results = BlockGf(name_block_generator = [ ( 'OpCorr',f(2) ) ], make_copies=False)

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
            for n,(proba,GM) in enumerate(self.global_moves):
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
        self.global_moves_mapping_list = []
        for n,(proba,GM) in enumerate(self.global_moves):
            F = extend_function_on_fundamentals(Map_GM_to_Fund_Ops(GM))
            m = {}
            for name,op in OPdict.items() :
                op_im = F(op)
                n1,n2 = myfind(op),myfind(op_im)
                m[n1] = n2
            name = "%s"%n
            self.global_moves_mapping_list.append((proba,m,name))
        #mpi.report ("global_moves_mapping_list", self.global_moves_mapping_list)

        # Now add the operator for F calculation if needed
        if self.use_f :
            Hloc_WithoutQuadratic = self.H_local.remove_quadratic()
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

        # Create the C_Cag_Ops array which describes the grouping of (C,Cdagger) operator
        # for the MonteCarlo moves : (a, alpha) block structure [ [ (C_name, Cdag_name)]]
        self.c_cdag_ops = [ [ (myfind(C(a,alpha)), myfind(Cdag(a,alpha))) for alpha in al ] for a,al in self.gf_struct]

        # Define G0_inv and correct it to have G0 to have perfect 1/omega behavior
        self.G0_inv = inverse(self.G0)
        Delta = self.G0_inv.copy()
        for n,d in Delta:
          d <<= A_Omega_Plus_B(self.G0_inv[n].tail[-1], self.G0_inv[n].tail[0])
          d -= self.G0_inv[n]
        for n,g in self.G0_inv:
          assert(g.N1==g.N2)
          identity=numpy.identity(g.N1)
          self.G0[n] <<= iOmega_n + g.tail[0]
          self.G0[n] -= Delta[n]
        self.G0_inv <<= self.G0
        self.G0.invert()

        # Construct the function in tau
        f = lambda g,L : GfImTime(indices = g.indices, beta = g.mesh.beta, n_points =L )
        self.Delta_tau = BlockGf(name_block_generator = [ (n,f(g,self.n_time_slices_delta) )   for n,g in self.G], make_copies=False, name='D')
        self.G_tau = BlockGf(name_block_generator = [ (n,f(g,self.n_time_slices_gtau) )    for n,g in self.G], make_copies=False, name='G')
        self.F_tau = BlockGf(name_block_generator = self.G_tau, make_copies=True, name='F')

        for (i,gt) in self.Delta_tau : gt.set_from_inverse_fourier(Delta[i])
        mpi.report("Inv Fourier done")
        if (self.legendre_accumulation):
            self.G_Legendre = BlockGf(name_block_generator = [ (n,GfLegendre(indices =g.indices, beta =g.mesh.beta, n_points =self.n_legendre) )   for n,g in self.G], make_copies=False, name='Gl')
        else:
            self.G_Legendre = BlockGf(name_block_generator = [ (n,GfLegendre(indices =[1], beta =g.mesh.beta, n_points =1) ) for n,g in self.G], make_copies=False, name='Gl') # G_Legendre must not be empty but is not needed in this case. So I make it as small as possible.

        # Starting the C++ code
        self.Sigma_old <<= self.Sigma

        # C++ solver
        solver_c(boost_object(self.__dict__),
                 new Hloc_c(nf, nb, boost_object(Oplist2), boost_object(QuantumNumberOperators), boost_object(SymList),
                            boost_object(self.quantum_number_selection), 0),
                 as_gf_block_imtime(self.G_tau),
                 as_gf_block_imtime(self.F_tau),
                 as_gf_block_imtime(self.Delta_tau),
                 as_gf_block_imtime(self.measured_time_correlators_results),
                 as_gf_block_legendre(self.G_Legendre)).solve()

        # Compute G on Matsubara axis possibly fitting the tail
        if self.legendre_accumulation:
          for s,g in self.G:
            identity=numpy.zeros([g.N1,g.N2],numpy.float)
            for i,m in enumerate (g.indicesL):
              for j,n in enumerate (g.indicesR):
                if m==n: identity[i,j]=1
            self.G_Legendre[s].enforce_discontinuity(identity)
            g <<= LegendreToMatsubara(self.G_Legendre[s])
        else:
          if (self.time_accumulation):
            for name, g in self.G_tau:
              identity=numpy.zeros([g.N1,g.N2],numpy.float)
              for i,m in enumerate (g.indicesL):
                for j,n in enumerate (g.indicesR):
                  if m==n: identity[i,j]=1
              g.tail.zero()
              g.tail[1] = identity
              self.G[name].set_from_fourier(g)

          # This is very sick... but what can we do???
          self.Sigma <<= self.G0_inv - inverse(self.G)
          self.fit_tails()
          self.G <<= inverse(self.G0_inv - self.Sigma)

        # Now find the self-energy
        self.Sigma <<= self.G0_inv - inverse(self.G)

        mpi.report("Solver %(name)s has ended."%self.__dict__)

        # for operator averages: if twice defined operator, rename output:
        for op1,op2 in self.twice_defined_ops.items():
            self.measured_operators_results[op1] = self.measured_operators_results[op2]
        for op1,op2 in self.twice_defined_ops.items():
            if op2 in self.measured_operators_results.keys(): del self.measured_operators_results[op2]

        if self.use_f :
            for (n,f) in self.F: f.set_from_fourier(self.F_tau[n])
            self.G2 = self.G0 + self.G0 * self.F
            self.Sigma2 = self.F * inverse(self.G2)

    #--------------------------------------------------

    def fit_tails(self):
        """ fits the tails of the noise monte carlo data """

        for n,sig in self.Sigma:

            known_coeff = numpy.zeros([sig.N1,sig.N2,1],numpy.float_)
            msh = [x.imag for x in sig.mesh]
            fit_start = msh[self.fit_start]
            fit_stop  = msh[self.fit_stop-1]

            sig.fit_tail(fixed_coef = known_coeff, order_max = 3, fit_start = fit_start, fit_stop = fit_stop)
