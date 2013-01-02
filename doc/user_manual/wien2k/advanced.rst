.. _advanced:

A more advanced example
=======================

Normally, one wants to adjust some more parameters in order to make the calculation more efficient. Here, we
will see a more advanced example, which is also suited for parallel execution. 
First, we load the necessary modules::

  from pytriqs.dft.sumk_lda import *
  from pytriqs.dft.converters.wien2k_converter import *
  from pytriqs.dft.solver_multiband import *
  from pytriqs.base.gf_local import *
  from pytriqs.base.archive import *

Then we define some parameters::

  LDAFilename='srvo3'
  U = 2.7
  J = 0.65
  Beta = 40
  Loops =  10                      # Number of DMFT sc-loops
  Mix = 1.0                        # Mixing factor of Sigma after solution of the AIM
  DeltaMix = 1.0                   # Mixing factor of Delta as input for the AIM
  DC_type = 1                      # DC type: 0 FLL, 1 Held, 2 AMF
  useBlocs = True                  # use bloc structure from LDA input
  useMatrix = False                # True: Slater parameters, False: Kanamori parameters U+2J, U, U-J
  use_spinflip = False             # use the full rotational invariant interaction?
  prec_mu = 0.0001
  QMCcycles = 20000
  Length_cycle = 200
  Warming_iterations = 2000

Most of these parameters are self-explaining. The first, `LDAFilename`, gives the filename of the input files. 
The next step, as described in the previous section, is to convert the input files::

  Converter = SumK_LDA_Wien2k_input(Filename=LDAFilename,repacking=True)
  Converter.convert_DMFT_input()
  mpi.barrier()

The command ``mpi.barrier()`` ensures that all nodes wait until the conversion of the input is finished on the master
node. After the conversion, we can check in the hdf5 archive, if previous runs are present, or if we have to start
from scratch::

  previous_runs = 0
  previous_present = False
  if mpi.is_master_node():
      ar = HDF_Archive(LDAFilename+'.h5','a')
      if 'iterations' in ar:
          previous_present = True
          previous_runs = ar['iterations']
      del ar
  previous_runs    = mpi.bcast(previous_runs)
  previous_present = mpi.bcast(previous_present)
  # if previous runs are present, no need for recalculating the bloc structure:
  calc_blocs = useBlocs and (not previous_present)

Now we can use all this information to initialise the :class:`SumK_LDA` class::

  SK=SumK_LDA(HDFfile=LDAFilename+'.h5',UseLDABlocs=calc_blocs)

If there was a previous run, we know already about the block structure, and therefore `UseLDABlocs` is set to `False`.
The next step is to initialise the Solver::

  Norb = SK.corr_shells[0][3]
  l = SK.corr_shells[0][2]
  S=Solver_MultiBand(Beta=Beta,U_interact=U,J_Hund=J,Norb=Norb,useMatrix=useMatrix, 
                     T=SK.T[0] ,GFStruct=SK.GFStruct_Solver[0],map=SK.map[0], 
                     l=l, deg_orbs=SK.deg_shells[0], use_spinflip=use_spinflip)

As we can see, many options of the solver are set by properties of the :class:`SumK_LDA` class, so we don't have
to set them manually. We now set the basic parameters of the QMC solver::

  S.N_Cycles  = QMCcycles
  S.Length_Cycle = Length_cycle
  S.N_Warmup_Cycles = Warming_iterations

If there are previous runs stored in the hdf5 archive, we can now load the self energy
of the last iteration::

  if (previous_present):
    if (mpi.is_master_node()):
        ar = HDF_Archive(LDAFilename+'.h5','a')
        S.Sigma <<= ar['SigmaF']
        del ar
    S.Sigma = mpi.bcast(S.Sigma)
    
The last command is the broadcasting of the self energy from the master node to the slave nodes. 
Now we can go to the definition of the self-consistency step. It consists again of the basic steps discussed in the 
previous section, with some additional refinement::

  for IterationNumber in range(1,Loops+1) :
     
        SK.symm_deg_GF(S.Sigma,orb=0)                           # symmetrise Sigma
        SK.put_Sigma(Sigmaimp = [ S.Sigma ])                    # put Sigma into the SumK class:

        Chemical_potential = SK.find_mu( precision = prec_mu )  # find the chemical potential
        S.G <<= SK.extract_Gloc()[0]                            # calculation of the local Green function
        mpi.report("Total charge of Gloc : %.6f"%S.G.total_density())

        if ((IterationNumber==1)and(previous_present==False)):
            # Init the DC term and the real part of Sigma, if no previous run was found:
            dm = S.G.density()
            SK.SetDoubleCounting( dm, U_interact = U, J_Hund = J, orb = 0, useDCformula = DC_type)
            S.Sigma <<= gf_init.Const(SK.dc_imp[0]['up'][0,0])
        
        # now calculate new G0:
        if (mpi.is_master_node()):
            # We can do a mixing of Delta in order to stabilize the DMFT iterations:
            S.G0 <<= S.Sigma + inverse(S.G)
            ar = HDF_Archive(LDAFilename+'.h5','a')
            if ((IterationNumber>1) or (previous_present)):
                mpi.report("Mixing input Delta with factor %s"%DeltaMix)
                Delta = (DeltaMix * S.G0.Delta()) + (1.0-DeltaMix) * ar['DeltaF']
                S.G0 <<= S.G0 + S.G0.Delta() - Delta
                
            ar['DeltaF'] = S.G0.Delta()
            S.G0 <<= inverse(S.G0)
            del ar
            
        S.G0 = mpi.bcast(S.G0)

        # Solve the impurity problem:
        S.Solve()

        # solution done, do the post-processing:
        mpi.report("Total charge of impurity problem : %.6f"%S.G.total_density())

        # Now mix Sigma and G with factor Mix, if wanted:
        if ((IterationNumber>1) or (previous_present)):
            if (mpi.is_master_node()):
                ar = HDF_Archive(LDAFilename+'.h5','a')
                mpi.report("Mixing Sigma and G with factor %s"%Mix)
                S.Sigma <<= Mix * S.Sigma + (1.0-Mix) * ar['SigmaF']
                S.G <<= Mix * S.G + (1.0-Mix) * ar['GF']
                del ar
            S.G = mpi.bcast(S.G)
            S.Sigma = mpi.bcast(S.Sigma)

        # Write the final Sigma and G to the hdf5 archive:
        if (mpi.is_master_node()):
            ar = HDF_Archive(LDAFilename+'.h5','a')
            ar['iterations'] = previous_runs + IterationNumber	
            ar['SigmaF'] = S.Sigma
            ar['GF'] = S.G
	    del ar

        # Now set new double counting:
        dm = S.G.density()
        SK.SetDoubleCounting( dm, U_interact = U, J_Hund = J, orb = 0, useDCformula = DC_type)
        
	#Save stuff:
        SK.save()

This is all we need for the LDA+DMFT calculation. At the end, all results are stored in the hdf5 output file.



