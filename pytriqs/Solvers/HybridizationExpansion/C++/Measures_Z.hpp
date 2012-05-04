
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

#ifndef TRIQS_CTHYB1_MEASURES_Z_H
#define TRIQS_CTHYB1_MEASURES_Z_H

#include "Python.h"

#include "triqs/mc_tools/mc_measure_set.hpp"
#include "boost/serialization/complex.hpp"

// Accumulate the sign : used as base of other classes
template<typename MCSignType>
class Measure_acc_sign { 
 uint64_t MeasureNumber;
 protected:
 typedef MCSignType mc_weight_type;
 MCSignType acc_sign;
 Measure_acc_sign() { acc_sign=0; MeasureNumber=0;}
 void accumulate(MCSignType signe) { MeasureNumber++; acc_sign += signe; }

 // after collect_results the acc_sign is reduced AND bcasted on all the nodes.
 void collect_results( boost::mpi::communicator const & c){
  uint64_t ntot=1;// to avoid the assert on the nodes...
  MCSignType acc_sign_tot=0;
  boost::mpi::reduce(c, acc_sign, acc_sign_tot, std::plus<MCSignType>(), 0);
  boost::mpi::reduce(c, MeasureNumber, ntot, std::plus<uint64_t>(), 0);
  assert (ntot>0);
  acc_sign = acc_sign_tot; // ok on master
  boost::mpi::broadcast(c, acc_sign,0);
  MeasureNumber =ntot;
  if (c.rank()==0) { 
   MCSignType sign_qmc = acc_sign / MCSignType(ntot); 
   std::cerr << "Average sign: " << abs(sign_qmc) << std::endl;
   if (abs(sign_qmc)<=1E-5) std::cerr << " MAJOR SIGN PB IN !"<<std::endl;
   if (abs(sign_qmc) < 0.01) std::cerr << "Very severe sign problem "<< sign_qmc<<std::endl;
  }  
 }
};

#endif

