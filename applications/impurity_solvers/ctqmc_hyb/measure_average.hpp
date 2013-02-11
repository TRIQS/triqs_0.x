
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

#ifndef TRIQS_CTHYB1_MEASURES_OpAv_H
#define TRIQS_CTHYB1_MEASURES_OpAv_H

#include "configuration.hpp"
#include "measure_z.hpp"

/**
   Measure the average of an operator
*/
class Measure_OpAv : public Measure_acc_sign<std::complex<double>> {
  Configuration & Config;
  double opAv;
  boost::python::dict opAv_res;
  typedef Measure_acc_sign<std::complex<double>> BaseType;  
public :   
  Measure_OpAv(std::string opName, Configuration & Config_, boost::python::dict opAv_results):
   BaseType(), Config(Config_), opAv(0.0), opAv_res(opAv_results) , name(opName){}

  const std::string name;

  void accumulate(std::complex<double> signe) {
   BaseType::accumulate(signe);
   const double s(real(signe)); // not elegant !
   double tau = Config.Beta/2.0;
   Config.DT.insertOneOperator(tau, Config.H[name]);
   double r = Config.DT.ratioNewTrace_OldTrace();
   Config.DT.undo_insertOneOperator();
   opAv += s*r;
  }

  void collect_results( boost::mpi::communicator const & c){
   BaseType::collect_results(c);
   double Z_qmc ( real(this->acc_sign));
   double op_average;
   boost::mpi::reduce(c, opAv, op_average, std::plus<double>(), 0);
   op_average /= Z_qmc;
   if (c.rank()==0) std::cout << "< " << name << " > = " << op_average << std::endl;
   opAv_res[name] = op_average; // master only
  }
};

#endif
