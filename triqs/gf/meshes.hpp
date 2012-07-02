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
#ifndef TRIQS_GF_LOCAL_MESHES_H
#define TRIQS_GF_LOCAL_MESHES_H

#include "domains.hpp"
#include "mesh_pt.hpp"

namespace triqs { namespace gf { namespace meshes { 

 struct matsubara_freq {

  typedef domains::matsubara_freq domain_type;
  typedef std::complex<double> gf_result_type;
  typedef long index_type;
  //typedef tqa::range range_type;
  //typedef arrays::range slice_arg_type;

  matsubara_freq (double Beta=1, statistic_enum s=Fermion, size_t N_max=1025): 
   _dom(Beta,s), n_max_(N_max), pi_over_beta(std::acos(-1)/Beta), sh(s==Fermion? 1:0) {}

  //domain_type domain() const { return _dom;}
  domain_type const & domain() const { return _dom;}
  size_t size() const{ return n_max_;}

  mesh_pt<matsubara_freq> operator[](index_type i) const { return make_mesh_pt(*this,i);}
 
  // REMOVE THIS EMBED !!
  domain_type::embedded_point_type embed(index_type n) const {return std::complex<double>(0,(2*n+sh)*pi_over_beta);}

  typedef  mesh_pt_generator<matsubara_freq> iterator;
  iterator begin() const { return iterator (*this);}
  iterator end()   const { return iterator (*this, true);}
 
  // move this out of the mesh 
  template<typename F> typename F::mv_type interpolate( F const & f, domain_type::point_type const & n) const {  
   return f((*this)[n]);
  } // here complex n -> -n and the tail : protection ...

  friend bool operator == (matsubara_freq M1, matsubara_freq M2) { return ((M1._dom == M2._dom) && (M1.n_max_ ==M2.n_max_) );}

  protected:
  domain_type _dom;
  size_t n_max_; 
  double pi_over_beta;
  long sh;
 };

 //--------------------------------------------------------

 class matsubara_time {

  public: 
  typedef domains::matsubara_time domain_type;
  typedef double gf_result_type;
  typedef long index_type;

  matsubara_time (double Beta=1, statistic_enum s=Fermion, size_t N_time_slices=1025, double delta=0):
   _dom(Beta,s), n_time_slices_(N_time_slices),_delta(delta), beta_over_l(Beta/N_time_slices){}

  domain_type const & domain() const { return _dom;}
  size_t size() const{ return n_time_slices_;}

  mesh_pt<matsubara_time> operator[](index_type i) const { return make_mesh_pt(*this,i);}

 // REMOVE THIS EMBED !!
  domain_type::embedded_point_type embed(index_type n) const {return _delta + n * beta_over_l;}

  typedef  mesh_pt_generator<matsubara_time> iterator;
  iterator begin() const { return iterator (*this);}
  iterator end()   const { return iterator (*this, true);}

  friend bool operator == (matsubara_time M1, matsubara_time M2) { return ((M1._dom == M2._dom) && (M1.n_time_slices_ ==M2.n_time_slices_) );}
 
  double delta() const { return _delta;}

  protected:
  domain_type _dom;
  size_t n_time_slices_;
  double _delta;
  double beta_over_l;
};

 struct real_freq {};

//--------------------------------------------------------

 struct BZ_mesh{
  typedef domains::BZ domain_type;
  typedef long index_type;
  typedef arrays::range slice_args_type;

  domain_type const & domain() const {return dom_;}
  mesh_pt<BZ_mesh> operator[](index_type n) const  { return make_mesh_pt(*this,n);}
  domain_type::point_type embed(index_type const & n) const {return 0;}

  template<typename F> typename F::mv_type interpolate( F const & f, domain_type::point_type const & k) const {   return f((*this)[ 0 ]); } //just to fill in the fct

  protected:
  domain_type dom_;
 };
}

}}

#endif




