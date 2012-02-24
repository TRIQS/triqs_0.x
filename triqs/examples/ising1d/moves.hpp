#ifndef moves_jadlsdadsjl
#define moves_jadlsdadsjl

#include <triqs/mc_tools/polymorphic_random_generator.hpp>
#include <vector>

// definition of a configuration
struct configuration {

  // N is the length of the chain, M the total magnetization
  // beta the inverse temperature, J, the magnetic field and the energy of the configuration
  int N, M;
  double beta, J, field, energy;

  // the chain of spins: true means "up", false means "down"
  std::vector<bool> chain;

  // constructor
  configuration(int N_, double beta_, double J_, double field_):
    N(N_), M(-N), beta(beta_), J(J_), field(field_), energy(-N*(J-field)), chain(N,false) {}

};


// flip a random spin
struct flip {

  typedef double mc_weight_type;
  configuration * config;
  triqs::mc_tools::polymorphic_random_generator &RNG;

  int site;
  double delta_energy;

  // constructor
  flip(configuration & config_, triqs::mc_tools::polymorphic_random_generator & RNG_) : config(&config_), RNG(RNG_) {}

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


// measure the magnetization
struct compute_m {

  configuration * config;
  int Z, M;

  compute_m(configuration & config_) : config(&config_), Z(0), M(0) {}

  // accumulate Z and magnetization
  void accumulate(int sign) {

    Z += sign;
    M += config->M;

  }

  // get final answer M / (Z*N)
  void collect_results(boost::mpi::communicator const &c) {

    std::cout << "Z = " << Z << std::endl;
    std::cout << "Magnetization: " << double(M)/(Z*config->N) << std::endl << std::endl;

  }

};

#endif
