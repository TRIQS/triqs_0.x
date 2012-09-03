.. highlight:: c

A simple Monte Carlo example
----------------------------

To be concrete, let's start with a full Monte Carlo code. We will consider
maybe the simplest problem ever: a single spin in a magnetic field :math:`h`
at a temperature :math:`1/\beta`. The Hamiltonian is simply:

.. math::

  \mathcal{H} = - h (n_\uparrow - n_\downarrow).

You want to compute the magnetization of this single spin. From statistical
mechanics it is clearly just

.. math::

  m = \frac{\exp(\beta h) - \exp(-\beta h)}{\exp(\beta h) + \exp(-\beta h)}


The C++ code for this problem
*****************************

Let's see how we can get this result from a Monte Carlo simulation. Here is
a code that would do the job. Note that we put everything in one file here,
but obviously you would usually want to cut this into pieces for clarity::

  #include <Python.h>
  #include <iostream>
  #include <boost/python.hpp>
  #include <triqs/mc_tools/mc_generic.hpp>

  // a configuration describing my system
  struct configuration {

    // the spin, the inverse temperature, the external field
    int spin; double beta, h;
    configuration(int spin_, double beta_, double h_) : spin(spin_), beta(beta_), h(h_) {}

  };

  // a single move: flip the spin
  struct flip {

    // the mc_weight_type is the type of the MC weight
    typedef double mc_weight_type;

    configuration & config;
    flip(configuration & config_) : config(config_) {}

    // The three methods below are expected for a move:
    // - Try() returns the Metropolis ratio
    // - Accept() updates the configuration and possibly corrects the sign
    // - Reject() does nothing in this case
    mc_weight_type Try() { return std::exp(-2*config.spin*config.h*config.beta); }
    mc_weight_type Accept() { config.spin *= -1; return 1.0; }
    void Reject() {}

  };

  // a measurement: the magnetization
  struct compute_m {

    configuration & config;
    int Z, M;

    compute_m(configuration & config_) : config(config_), Z(0), M(0) {}

    // These two methods are expected for a measurement:
    // - accumulate(s) is called at every measure
    // - collect_results(c) is called at the very end
    void accumulate(int sign) { Z += sign; M += config.spin; }
    void collect_results(boost::mpi::communicator const &c) {
      std::cout << "Magnetization: " << double(M)/Z << std::endl << std::endl; }

  };

  // The main code
  int main() {

    boost::mpi::communicator c;

    // greeting
    std::cout << "Isolated spin" << std::endl << std::endl;

    // Initialize the python. Otherwise no boost::python will work
    Py_Initialize();

    // Prepare the MC parameters in a python dict
    boost::python::dict d;
    d["N_Cycles"] = 5000000;
    d["Length_Cycle"] = 10;
    d["N_Warmup_Cycles"] = 10000;
    d["Random_Seed"] = 374982;
    d["Verbosity"] = 1;

    // Construct a Monte Carlo loop
    triqs::mc_tools::mc_generic<double> SpinMC(d, 0);

    // parameters of the model
    double beta = 0.3;
    double field = 0.5;

    // construct configuration
    configuration config(-1, beta, field);

    // add moves and measures
    SpinMC.add_move(new flip(config), 1.0, "spin flip");
    SpinMC.add_measure(new compute_m(config));

    // Run and collect results
    SpinMC.run(triqs::mc_tools::clock_callback(-1));
    SpinMC.collect_results(c);

    // Finalize everything
    Py_Finalize();
    return 0;

  }

Let's go through the different parts of this code


Setting the Monte Carlo parameters
**********************************

Let's start with the ``main()``. Before constructing the Monte Carlo
class we need to gather the parameters in a dictionary that will
be passed as a parameter at construction. This is done with these
lines::

    boost::python::dict d;
    d["N_Cycles"] = 5000000;
    d["Length_Cycle"] = 10;
    d["N_Warmup_Cycles"] = 10000;
    d["Random_Seed"] = 374982;
    d["Verbosity"] = 1;

The keys ``N_Cycles``, ``Length_Cycle`` and ``N_Warmup_Cycles`` determine the
length of the Monte Carlo cycles, the number of measurements and the warmup
length. This has been detailed earlier. The parameter ``Random_Seed`` sets the
seed for the random number generator.  Finally ``Verbosity`` sets the verbosity
level. There's some output if it is 1 and essentially no output if it is 0.
All these parameters are mandatory.

Constructing the Monte Carlo simulation
***************************************

The Monte Carlo simulation is constructed with::

    triqs::mc_tools::mc_generic<double> SpinMC(d, 0);

Note that you need to include the header ``<triqs/mc_tools/mc_generic.hpp>``
in order to access the ``mc_generic`` class. The simulation is
constructed from the parameter dictionary and a 0.

Moves and measures
******************

At this stage the basic structure of the Monte Carlo is in ``SpinMC``.  But we
now need to tell it what moves must be tried and what measures must be made.
This is done with::

    SpinMC.add_move(new flip(config), 1.0, "spin flip");
    SpinMC.add_measure(new compute_m(config));

The method ``add_move`` expects a pointer to a move, a number we'll explain
later and a name. The measure expects a pointer to a measure. As you can
see we add a "flip" move and a measurement of the magnetization.




