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

  struct tail{
   typedef std::complex<double> gf_result_type;
   static const bool has_tail = false;
  };

  struct matsubara_freq {
   typedef std::complex<double> gf_result_type;
   double beta;
   statistic_enum statistic;
   matsubara_freq (double Beta, statistic_enum s = Fermion): beta(Beta), statistic(s){}
   static const bool has_tail = true;
  };

  struct matsubara_time {
   typedef double gf_result_type;
   double beta;
   statistic_enum statistic;
   matsubara_time (double Beta, statistic_enum s = Fermion): beta(Beta), statistic(s){}
   static const bool has_tail = true;
  };

  struct matsubara_legendre{};
  struct real_freq {};
  struct real_time {};
 }

//#define TRIQS_LOCAL_GF_DOMAIN_LIST (matsubara_freq)(matsubara_time)(matsubara_legendre)(real_freq)(real_time)

 namespace meshes { 

  class tail{
   int omin,omax;
 
   public:
   typedef domains::tail domain_type;
   typedef std::complex<double> gf_result_type;
   typedef size_t index_type;
   typedef tqa::range range_type;
   
   tail(int OrderMin=-1, int OrderMax = 5) : omin(OrderMin), omax(OrderMax) {}
   
   static const bool has_tail = false;
   static const bool mesh_tail = false;
   
   size_t size() const{ int r = omax - omin +1; assert(r>=0); return r;}
   int order_min() const {return omin;}
   int order_max() const {return omax;}
  };

  //--------------------------------------------------------

  struct matsubara_freq {
 
   typedef domains::matsubara_freq domain_type;
   typedef std::complex<double> gf_result_type;
   typedef size_t index_type;
   typedef tqa::range range_type;

   matsubara_freq (double Beta=1, statistic_enum s=Fermion, size_t N_max=1025, size_t tail_expansion_order=5 ): 
    mesh_tail(tail_expansion_order), _dom(Beta,s), n_max_(N_max) {}

   tail mesh_tail;
   
   domain_type const & domain() const { return _dom;}
   size_t size() const{ return n_max_;}

   protected:
   domain_type _dom;
   size_t n_max_; 
  };

  //--------------------------------------------------------
/*
  class matsubara_time {
   size_t n_time_slices_;
 
   public: 
   typedef domains::matsubara_time domain_type;
   typedef double gf_result_type;
   typedef size_t index_type;
   
   matsubara_time (size_t N_time_slices, Statistic s): n_time_slices_( N_time_slices), statistic(s) {}
   
   //size_t n_time_slices() const { return n_time_slices_;}
   
   Statistic statistic;
   
   static const bool has_tail = true;
   tail mesh_tail;
   
   size_t size() const{ return n_time_slices_;}
  };

  */
 }

}}

#endif




