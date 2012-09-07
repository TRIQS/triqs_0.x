#include <Python.h>
#include <iostream>
#include <boost/python.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include "moves.hpp"

int main() {

  boost::mpi::communicator c;

  // greeting
  std::cout << "Ising 1D" << std::endl << std::endl;

  // Initialize the python. Otherwise no boost::python will work
  Py_Initialize();
  
  // Prepare the MC parameters
  boost::python::dict d;
  d["N_Cycles"] = 500000;
  d["Length_Cycle"] = 50;
  d["N_Warmup_Cycles"] = 10000;
  d["Random_Seed"] = 374982;
  d["Verbosity"] = 1;

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<double> IsingMC(d, 0);

  // parameters of the model
  int length = 100;
  double J = -1.0;
  double field = 0.5;
  double beta = 0.3;

  // construct configuration
  configuration config(length, beta, J, field);

  // add moves and measures
  IsingMC.add_move(new flip(config, IsingMC.RandomGenerator), 1.0, "spin flip");
  IsingMC.add_measure(new compute_m(config));

  // Run and collect results
  IsingMC.run(triqs::utility::clock_callback(-1));
  IsingMC.collect_results(c);

  // Finalize everything
  Py_Finalize();
  return 0;

}
