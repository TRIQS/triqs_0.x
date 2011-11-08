The Measure concept
===================

* **Purpose**  : Defines the basic MonteCarlo measure.
* **Definition** : 

  ==========================================================================  ============================================================
  Elements                                                                    Comment
  ==========================================================================  ============================================================
  * void accumulate(std::complex<double> sign)                                - Accumulation with the sign
  * void collect_results ( boost::mpi::communicator const & c)                - Collects the results over the communicator, and finalize
                                                                                the calculation (compute average, error). 
  ==========================================================================  ============================================================



