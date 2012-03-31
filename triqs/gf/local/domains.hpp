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

 enum Statistic {Boson,Fermion};

 namespace domains {
  struct matsubara_freq{};
  struct matsubara_time{};
  struct matsubara_legendre{};
  struct real_freq {};
  struct real_time {};
  struct infty{};
 }

 namespace meshes { 

  struct tail{
   typedef std::complex<double> gf_result_type;
   static const bool has_tail = false;
   typedef size_t index_type;
   static const bool mesh_tail = false;
   const size_t order;
   size_t len() const{ return order;}
   tail(size_t order_=5) : order(order_) {}
  };

  struct matsubara_freq : domains::matsubara_freq {
   const size_t n_max; 
   Statistic statistic;
   static const bool has_tail = true;
   tail mesh_tail;
   typedef std::complex<double> gf_result_type;
   typedef size_t index_type;
   size_t len() const{ return n_max;}
   matsubara_freq (Statistic s=Fermion, size_t n_max_=1025, size_t tail_expansion_order=5 ): 
    n_max( n_max_), statistic(s), mesh_tail(tail_expansion_order) {}
  };

  struct matsubara_time : domains::matsubara_freq {
   const size_t n_time_slices; 
   Statistic statistic;
   static const bool has_tail = true;
   typedef double gf_result_type;
   typedef size_t index_type;
   size_t len() const{ return n_time_slices;}
   matsubara_time (size_t n_time_slices_, Statistic s): n_time_slices( n_time_slices_), statistic(s) {}
  };

  
 }

}}

#endif




