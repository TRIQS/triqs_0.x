/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_LOCAL_DOMAIN_H
#define TRIQS_GF_LOCAL_DOMAIN_H

namespace triqs { namespace gf {

 namespace tqa= triqs::arrays;

 enum statistic_enum {Boson,Fermion};

 namespace domains {

  struct infty{}; // the point at infinity
  
  struct matsubara_freq {
   typedef long                 point_type;
   typedef std::complex<double> embedded_point_type;
   typedef std::complex<double> gf_result_type;
   double beta;
   statistic_enum statistic;
   matsubara_freq (double Beta, statistic_enum s = Fermion): beta(Beta), statistic(s){}
   static const bool has_tail = true;
   bool operator == (matsubara_freq D) { return ((std::abs(beta - D.beta)<1.e-15) && (statistic == D.statistic));}
  };

  struct matsubara_time {
   typedef double point_type;
   typedef double embedded_point_type;
   typedef double gf_result_type;
   double beta;
   statistic_enum statistic;
   matsubara_time (double Beta, statistic_enum s = Fermion): beta(Beta), statistic(s){}
   //static const bool has_tail = true;
  };

  struct matsubara_legendre{};
  struct real_freq {};
  struct real_time {};
 }

 //#define TRIQS_LOCAL_GF_DOMAIN_LIST (matsubara_freq)(matsubara_time)(matsubara_legendre)(real_freq)(real_time)

}}

#endif




