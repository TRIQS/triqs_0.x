
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

#ifndef TRIQS_CTHYB1_MEASURES_OpCorr_H
#define TRIQS_CTHYB1_MEASURES_OpCorr_H

#include <triqs/gf/imtime.hpp>
#include "measure_z.hpp"
#include "configuration.hpp"
#include "gf_binner.hpp"

/**
   Measure the time correlator of an operator: <O(tau)O>
*/
class Measure_OpCorr : public Measure_acc_sign<double> {
  const std::string opName;
  Configuration & Config;
  const int N_timeslices;
  const double deltatau;
  triqs::gf::gf_view<triqs::gf::imtime> Op_res;
  gf_binner<triqs::gf::gf_view<triqs::gf::imtime>> Op_res_bin;
  typedef Measure_acc_sign<double> BaseType;
public :   
  Measure_OpCorr(std::string MeasureName_, std::string opName_, Configuration & Config_, triqs::gf::gf_view<triqs::gf::imtime> &Op_res_, int N_timeslices_):
    BaseType(),  opName(opName_), Config(Config_),  N_timeslices(N_timeslices_),deltatau(Config.Beta/N_timeslices_),Op_res(Op_res_),Op_res_bin(Op_res), 
     name(opName)  {}
 
  const std::string name;
  
  void accumulate(std::complex<double> signe);  

  void collect_results( boost::mpi::communicator const & c){
  BaseType::collect_results(c);
   mc_weight_type Z_qmc ( this->acc_sign);

   auto res = triqs::make_clone(Op_res);
   auto g_loc = triqs::make_clone(Op_res);
   boost::mpi::reduce(c, g_loc, res, std::plus<triqs::gf::gf<triqs::gf::imtime>>(),0);
   boost::mpi::broadcast(c,res,0);
   Op_res = res / (Z_qmc * deltatau);

  }
  
};


#endif
