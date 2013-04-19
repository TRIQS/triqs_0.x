
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef TRIQS_CTHYB1_MEASURES_G_H
#define TRIQS_CTHYB1_MEASURES_G_H

#include <triqs/gf/imtime.hpp>
#include "configuration.hpp"
#include "measure_z.hpp"
#include "gf_binner.hpp"

/**
   Measure the Green function (one bloc only) in time.
*/
class Measure_G_tau : public Measure_acc_sign<std::complex<double>> {
  typedef Measure_acc_sign<std::complex<double>> BaseType;
  const std::string name;  
  const Configuration & Config;
  triqs::gf::gf_view<triqs::gf::imtime> & G_tau;
  gf_binner<triqs::gf::gf_view<triqs::gf::imtime>> G_tau_bin;
  const int a_level;
 public :   

  Measure_G_tau(const Configuration & Config_,int a, triqs::gf::gf_view<triqs::gf::imtime> & Gtau_):
   BaseType(), name(to_string("G(tau)",a)), Config(Config_), G_tau(Gtau_), G_tau_bin(G_tau), a_level(a) { }

  void accumulate(std::complex<double> signe)  { 
   BaseType::accumulate(signe);
   const double s(real(signe)); // not elegant !
   for (Configuration::DET_TYPE::C_Cdagger_M_iterator p(*Config.dets[a_level]); !p.atEnd(); ++p) 
    G_tau_bin(Config.info[p.C()->Op->Number].alpha, p.C()->tau, 
      Config.info[p.Cdagger()->Op->Number].alpha, p.Cdagger()->tau, 
      s* p.M());
  }

  void collect_results( boost::mpi::communicator const & c){

   BaseType::collect_results(c);
   mc_weight_type Z_qmc ( this->acc_sign);

   auto res = triqs::make_clone(G_tau);
   auto loc_g = triqs::make_clone(G_tau);
   boost::mpi::reduce(c, loc_g, res, std::plus<triqs::gf::gf<triqs::gf::imtime>>(),0);
   boost::mpi::broadcast(c,res,0);
   G_tau = res / (- real(Z_qmc) * Config.Beta * G_tau.mesh().delta());

  }

};

#endif
