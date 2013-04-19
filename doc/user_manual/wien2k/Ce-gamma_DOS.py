from pytriqs.applications.dft.sumk_lda_tools import *
from pytriqs.applications.dft.converters.wien2k_converter import *
from pytriqs.applications.impurity_solvers.hubbard_I.solver import Solver

# Creates the data directory, cd into it:
#Prepare_Run_Directory(DirectoryName = "Ce-Gamma") 
LDAFilename = 'Ce-gamma'
Beta =  40
Uint = 6.00
JHund = 0.70
DC_type = 0                      # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
load_previous = True              # load previous results
useBlocs = False                 # use bloc structure from LDA input
useMatrix = True                 # use the U matrix calculated from Slater coefficients instead of (U+2J, U, U-J)
ommin=-4.0
ommax=6.0
N_om=2001
broadening = 0.02

HDFfilename = LDAFilename+'.h5'

# Convert DMFT input:
# Can be commented after the first run
Converter = Wien2kConverter(filename=LDAFilename,repacking=True)
Converter.convert_dmft_input()
Converter.convert_par_proj_input()

#check if there are previous runs:
previous_runs = 0
previous_present = False

if mpi.is_master_node():
    ar = HDFArchive(HDFfilename,'a')
    if 'iterations' in ar:
        previous_present = True
        previous_runs = ar['iterations']
    else: 
        previous_runs = 0
        previous_present = False
    del ar

mpi.barrier()
previous_runs    = mpi.bcast(previous_runs)
previous_present = mpi.bcast(previous_present)

# if previous runs are present, no need for recalculating the bloc structure
# It has to be commented, if you run this script for the first time, starting
# from a converted h5 archive.

# Init the SumK class
SK = SumkLDATools(hdf_file=LDAFilename+'.h5',use_lda_blocks=False)


if (mpi.is_master_node()):
    print 'DC after reading SK: ',SK.dc_imp[SK.invshellmap[0]]

N = SK.corr_shells[0][3]
l = SK.corr_shells[0][2]

# Init the Solver:
S = Solver(Beta = Beta, Uint = Uint, JHund = JHund, l = l)
S.Nmoments=10

# set atomic levels:
eal = SK.eff_atomic_levels()[0]
S.set_atomic_levels( eal = eal )
S.GF_realomega(ommin=ommin, ommax = ommax, N_om=N_om)
S.Sigma.save('S.Sigma')
SK.put_Sigma(Sigmaimp = [S.Sigma])
SK.dos_partial(broadening=broadening)
