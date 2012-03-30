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

#ifndef TRIQS_GF_LOCAL_GF_H
#define TRIQS_GF_LOCAL_GF_H 

#include <triqs/lazy/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays;
 namespace tql= triqs::lazy;
 using tqa::range;

 template<typename DomainType, bool IsView> class gf {

 public:

  typedef double value_element_type;
  typedef DomainType domain_type;
  typedef tqa::matrix_view<element_value_type, Option::Fortran>       result_type;
  typedef tqa::matrix_view<const element_value_type, Option::Fortran> const_result_type;

  typedef void has_view_type_tag; // Do this or specialize the tempalte ????
  typedef gf_view<DomainType> view_type;
  typedef gf<DomainType>      non_view_type;

 protected:
  typedef tqa::array      <value_element_type,3>              data_type;
  typedef tqa::array_view <value_element_type,3>              data_view_type;
  typename mpl::if_<IsView, data_view_type, data_type>::type  data;

  std::vector<std::vector<string> > indices;
  mesh m_;

 public:
  
  gf (data_view_type const & data_, std::vector<string> const & indicesL, std::vector<string> const & indicesR) : 
   data(data_) {
    indices.push_back(indicesL); indices.push_back(indicesR);
   }

  gf (size_t N1, size_t N2, GridType const & g_, std::vector<string> const & indicesL, std::vector<string> const & indicesR) : 
   grid(g_), data(N1,N2,g_.len()), dom(d_) {
    indices.push_back(indicesL); indices.push_back(indicesR);
   }


  template<typename T>
   typename boost::disable_if< tql::is_lazy<T>, result_type>::type 
   operator() (T const & x) { return data(range(),range(),x);}
  
  template<typename T>
   typename boost::disable_if< tql::is_lazy<T>, const_result_type>::type 
   operator() (T const & x) const { return data(range(),range(),x);}

  TRIQS_LAZY_ADD_LAZY_CALL_WITH_COPY(1,gf);

  // pb : view of the indices !!
  template<typename A, typename B> view_type slice(A a, B b) { return gf_view( data( range(a), range(b), range())); }
  template<typename A, typename B> const view_type slice(A a, B b) const { return gf_view( data( range(a), range(b), range())); }

  view_type view() { return gf_view(data);}
  const view_type view() const { return gf_view(data);}

  template<typename RHS> // specialize for various RHS ( fourier_impl, other gf, etc....)
  gf & operator = (RHS const & rhs) { triqs::gf::local::assignment_impl<gf,RHS>::invoke(*this,rhs);} 

  // lazy_assignable
  template<typename F>
  void set_from_function(F f) { const size_t Nmax = data.shape()[2]; for (size_t u=0; u<Nmax; ++u) data(range(),range(),u) = F(u);}

  /// Save the Green function in i omega_n (as 2 columns).
  void save(string file,  bool accumulate=false) const;

  /// Load the GF
  void load(string file);

  /// HDF5 saving ....

};

// gf.impl.hpp

// implementation of expression template as function of omega...
// keep the domain

// 
// integrate + enable ??
template< typename GfType> struct integrate_impl; 
template< typename GfType> typename GfType::result_type integrate(GfType const & x) { return integrate_impl<GfType>::invoke(x);}
// + lazy overload
// 
// density + enable ??
template< typename GfType> struct density_impl; 
template< typename GfType> typename GfType::result_type density(GfType const & x) { return density_impl<GfType>::invoke(x);}
//+ lazy overload.

// 1 fichier de specialization par type de GF
//

// rewrite the immutable_concept : call, slice, domain

// implementation of save/load

}}

namespace lazy { 

  // overload of the free functions
  


}

}
#endif




