.. highlight:: c

.. _ising_solution:

The Ising chain in a magnetic field
-----------------------------------

Here is the a simple Monte-Carlo for a one-dimensional Ising chain.  The
problem is described in detail in this section about :ref:`the Ising model
<isingex>`.

The configuration
*****************

We start by defining a configuration class on which the move and measure
classes will act. We write this class in a file :file:`configuration.hpp`::

    #ifndef configuration_hpp
    #define configuration_hpp

    // The configuration of the system
    struct configuration {

      // N is the length of the chain, M the total magnetization
      // beta the inverse temperature, J the coupling , field the magnetic field and energy the energy of the configuration
      // field the magnetic field and energy the energy of the configuration
      int N, M;
      double beta, J, field, energy;

      // the chain of spins: true means "up", false means "down"
      std::vector<bool> chain;

      // constructor
      configuration(int N_, double beta_, double J_, double field_):
        N(N_), M(-N), beta(beta_), J(J_), field(field_), energy(-N*(J-field)), chain(N,false) {}

    };

    #endif

The move
********

The 
The move class should have three methods: `Try()`, `Accept()` and `Reject()`::

    #ifndef moves_hpp
    #define moves_hpp

    #include <triqs/mc_tools/random_generator.hpp>
    #include <vector>
    #include "configuration.hpp"

    // A move flipping a random spin
    struct flip {

      configuration * config;
      triqs::mc_tools::random_generator &RNG;

      int site;
      double delta_energy;

      // constructor
      flip(configuration & config_, triqs::mc_tools::random_generator & RNG_) :
         config(&config_), RNG(RNG_) {}

      double Try() {
        // pick a random site
        site = RNG(config->N);

        // find the neighbors with periodicity
        int left = (site==0 ? config->N-1 : site-1);
        int right = (site==config->N-1 ? 0 : site+1);

        // compute energy difference from field
        delta_energy = (config->chain[site] ? 2 : -2) * config->field;

        // compute energy difference from J
        if(config->chain[left] == config->chain[right]) {
          delta_energy += (config->chain[left] == config->chain[site] ? 4 : -4) * config->J;
        }

        // return Metroplis ratio
        return std::exp(-config->beta * delta_energy);
      }

      // if move accepted just flip site and update energy and magnetization
      double Accept() {
        config->M += (config->chain[site] ? -2 : 2);
        config->chain[site] = !config->chain[site];
        config->energy += delta_energy;

        return 1.0;
      }

      // nothing to do if the move is rejected
      void Reject() {}
    };


    #endif


Measure
*******

The measure class has two methods, `accumulate` and `collect_results`::


      #ifndef MEASURES_HPP
      #define MEASURES_HPP

      #include "configuration.hpp"
      #include <triqs/arrays/h5/array_stack.hpp>

      namespace tqa=triqs::arrays;

      ////////////////// measure the magnetization ////////////////// 
      struct compute_m {

        configuration * config;
        int Z, M;
        tqa::h5::H5File outfile;
        tqa::h5::array_stack<tqa::array<double,1> > S;

        template<class T1>  //convert to string
         static std::string filename(T1 x1) { std::stringstream f; f<<x1; return f.str(); }
             
        compute_m(configuration & config_) :
          config(&config_), Z(0), M(0),
          outfile(filename("M_stack.h5"), H5F_ACC_TRUNC ),
          S(outfile, "M",tqa::make_shape(1), 10)
          {    std::cout<< "Measure ok"<<std::endl;    }

        // accumulate Z and magnetization
        void accumulate(int sign) {
          Z += sign;
          M += config->M;
          S << double(M)/(Z*config->N);
        }

        // get final answer M / (Z*N)
        void collect_results(boost::mpi::communicator const &c) {
          std::cout << "Z = " << Z << std::endl;
          std::cout << "Magnetization: " << double(M)/(Z*config->N) << std::endl << std::endl;
        }  
      };

      #endif



Main program
************

The Monte-Carlo itself can now be written::

    #include <Python.h>
    #include <iostream>
    #include <boost/python.hpp>
    #include <triqs/mc_tools/mc_generic.hpp>
    #include <triqs/utility/callbacks.hpp>

    #include "moves.hpp"
    #include "configuration.hpp"
    #include "measures.hpp"

    int main() {

      boost::mpi::communicator c;

      // Initialize the python. Otherwise no boost::python will work
      Py_Initialize();
      
      // Prepare the MC parameters
      boost::python::dict d;
      d["N_Cycles"] = 500000;
      d["Length_Cycle"] = 50;
      d["N_Warmup_Cycles"] = 100000;
      d["Random_Seed"] = 374982;
      d["Verbosity"] = 1;

      // Construct a Monte Carlo loop
      triqs::mc_tools::mc_generic<double> IsingMC(d, 0);

      // parameters of the model
      int length = 100;
      double J = -1.0;
      double field = 0;
      double beta = 0.3;

      // construct configuration
      configuration config(length, beta, J, field);

      // add moves and measures
      IsingMC.add_move(new flip(config, IsingMC.RandomGenerator), 1.0, "spin flip");
      std::cout << "Add measure"<<std::endl;
      IsingMC.add_measure(new compute_m(config));
      std::cout << "Run"<<std::endl;

      // Run and collect results
      IsingMC.start(1.0, triqs::utility::clock_callback(-1));
      IsingMC.collect_results(c);

      // Finalize everything
      Py_Finalize();
      return 0;
    }


This yields::

    Add measure
    Measure ok
    Run
    1%; 2%; 3%; 4%; 5%; 6%; 7%; 8%; 9%; 10%; 11%; 12%; 13%; 14%; 15%; 16%; 17%; 18%; 19%; 20%; 21%;
    22%; 23%; 24%; 25%; 26%; 27%; 28%; 29%; 30%; 31%; 32%; 33%; 34%; 35%; 36%; 37%; 38%; 39%; 40%; 41%;
    42%; 43%; 44%; 45%; 46%; 47%; 48%; 49%; 50%; 51%; 52%; 53%; 54%; 55%; 56%; 57%; 58%; 59%; 60%; 61%;
    62%; 63%; 64%; 65%; 66%; 67%; 68%; 69%; 70%; 71%; 72%; 73%; 74%; 75%; 76%; 77%; 78%; 79%; 80%; 81%;
    82%; 83%; 84%; 85%; 86%; 87%; 88%; 89%; 90%; 91%; 92%; 93%; 94%; 95%; 96%; 97%; 98%; 99%; Z = 500000
    Magnetization: -0.00025888

