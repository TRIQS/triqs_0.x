import pytriqs.utility.mpi as mpi
import copy
from ctqmc_solver import Solver

class SolverOld(Solver):

    def __init__(self, Beta, GFstruct, N_Matsubara_Frequencies=1025, **param):

        Solver.__init__(self, beta=Beta, gf_struct=GFstruct, n_matsubara=N_Matsubara_Frequencies)
        self.params = param
        self.gen_keys = copy.deepcopy(self.__dict__)

        msg = """
**********************************************************************************
 Warning: You are using the old constructor for the solver. Beware that this will
 be deprecated in future versions. Please check the documentation.
**********************************************************************************
"""
        mpi.report(msg)

    def Solve(self):

        self.params.update(self.__dict__)
        self.update_params(self.params)
        for i in self.gen_keys: self.params.pop(i)
        self.params.pop("gen_keys")
        self.solve(**self.params)

    # This is for backward compatibility
    def update_params(self, d):

      allparams = [
        ('QMC_N_cycles_MAX', 'n_cycles'),
        ('NCycles', 'n_cycles'),
        ('Hloc', 'H_local'),
        ('QuantumNumbers', 'quantum_numbers'),
        ('Length_One_QMC_Cycle', 'length_cycle'),
        ('Number_Warming_Iteration', 'n_warmup_cycles'),
        ('Number_Frequencies_Accumulated', 'fit_stop'),
        ('Global_Move', 'global_moves'),
        ('UseSegmentPicture', 'use_segment_picture'),
        ('Proba_Move_Insert_Remove_Kink', 'prob_insert_remove'),
        ('Proba_Move_Move_Kink', 'prob_move'),
        ('OperatorsToAverage', 'measured_operators'),
        ('OpCorrToAverage', 'measured_time_correlators'),
        ('RecordStatisticConfigurations', 'record_stat'),
        ('N_Cycles', 'n_cycles'),
        ('H_Local', 'H_local'),
        ('Quantum_Numbers', 'quantum_numbers'),
        ('Random_Generator_Name', 'random_name'),
        ('Random_Seed', 'random_seed'),
        ('Length_Cycle', 'length_cycle'),
        ('Legendre_Accumulation', 'legendre_accumulation'),
        ('N_Legendre_Coeffs', 'n_legendre'),
        ('Time_Accumulation', 'time_accumulation'),
        ('N_Warmup_Cycles', 'n_warmup_cycles'),
        ('N_Frequencies_Accumulated', 'fit_stop'),
        ('Fitting_Frequency_Start', 'fit_start'),
        ('N_Time_Slices_Delta', 'n_time_slices_delta'),
        ('N_Time_Slices_Gtau', 'n_time_slices_gtau'),
        ('Nmax_Matrix', 'n_max_matrix'),
        ('Global_Moves', 'global_moves'),
        ('Use_Segment_Picture', 'use_segment_picture'),
        ('Proba_Insert_Remove', 'prob_insert_remove'),
        ('Proba_Move', 'prob_insert_move'),
        ('Measured_Operators', 'measured_operators'),
        ('Measured_Time_Correlators', 'measured_time_correlators'),
        ('Record_Statistics_Configurations', 'record_stat'),
        ('Use_F', 'use_f'),
        ('Quantum_Numbers_Selection', 'quantum_number_selection'),
        ('MAX_TIME', 'max_time')
      ]

      for (old, new) in allparams:
        if old in d:
          val = d.pop(old)
          d.update({new:val})
