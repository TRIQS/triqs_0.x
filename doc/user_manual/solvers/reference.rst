.. _ctqmc_ref:

The Solver class
=================

.. autoclass:: pytriqs.solvers.ctqmc_hyb.Solver
  :members: 

 
Compulsory entries
--------------------

These are the compulsory entries to put in the constructor of the solver:

  =================  ==========  ============================================================
  Key                Type        Documentation
  =================  ==========  ============================================================
  N_Cycles           int         Maximum number of Monte Carlo cycles
  H_Local            instance    Local Hamiltonian
  Quantum_Numbers    dict        Quantum Numbers
  =================  ==========  ============================================================

Optional entries
-------------------

  These instead are optional and come with a default value if they are not set:

  ==================================  ==================================  ==========  ============================================================
  Key                                 Default                             Type        Documentation
  ==================================  ==================================  ==========  ============================================================
  Length_Cycle                        200                                 int         Length of one QMC cycle
  N_Warmup_Cycles                     1000                                int         Number of warming cycles
  N_Legendre_Coeffs                   30                                  int         Number of Legendre coefficients that are used in practice
  Use_Segment_Picture                 False                               bool        Guarantee is made that the C^+ C alternate in each config
  Random_Generator_Name               ""                                  str         Name of the random number generator
  Random_Seed                         34788                               int         Seed for the random generator
  Global_Moves                        []                                  list        A list of global moves
  Measured_Operators                  {}                                  dict        A dict of operators that will be averaged
  Measured_Time_Correlators           {}                                  dict        A dict of operators, whose time correlations are to be measured
  Proba_Move                          1.0                                 float       Probability to move operators
  Proba_Insert_Remove                 1.0                                 float       Probability to insert/remove operators
  N_Time_Slices_Gtau                  10000                               int         Number of times slices in G_tau
  N_Time_Slices_Delta                 10000                               int         Number of times slices in Delta
  Legendre_Accumulation               True                                bool        Do we accumulate in legendre?
  Time_Accumulation                   False                               bool        Do we accumulate in imaginary-time?
  N_Frequencies_Accumulated           100                                 int         Number of frequencies to accumulate
  Fitting_Frequency_Start             50                                  int         Frequency at which the fit starts
  Nmax_Matrix                         100                                 int         Initial size of the determinant matrices
  Eta                                 0.0                                 float       (Expert only) Value of eta, the minimum value of the det
  Record_Statistics_Configurations    False                               bool        (Expert only) Get the kink length statistics
  Use_F                               False                               bool        (Expert only) Compute F
  Keep_Full_MC_Series                 False                               bool        (Expert only) Store the Green's function for later analysis
  Quantum_Numbers_Selection           <function <lambda> at 0x29de9b0>    function    (Prototype) A function to select quantum numbers
  ==================================  ==================================  ==========  ============================================================
