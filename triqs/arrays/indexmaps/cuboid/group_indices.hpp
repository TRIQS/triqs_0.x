/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_GROUP_INDICES_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_GROUP_INDICES_H
#include "boost/mpl/max_element.hpp"
#include <boost/mpl/size.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <triqs/arrays/array.hpp>

namespace triqs { namespace arrays {

 // MIndexList is an MPL vector of m_index describing the grouping of indices
 template<typename A, typename MIndexList> struct group_indices_impl { 

  typedef typename A::indexmap_type::index_order_type A_index_order_type;
  static const int new_dim = boost::mpl::size<MIndexList>::value;

  // First from a group_of_indices, find their position in memory, checking that indices are contiguous
  struct one_index_to_memory_pos { 
   template <typename T> struct apply {
    typedef boost::mpl::int_<A_index_order_type::template index_to_memory_rank<T::value>::value> type; 
   };
  };

  struct one_index_group_to_memory_pos { 
   template<typename V1> struct apply {
    typedef typename V1::type V;
    typedef typename boost::mpl::transform <V, one_index_to_memory_pos >::type type;
    static const int mi = boost::mpl::deref< typename boost::mpl::min_element<type>::type  >::type::value;	
    static const int ma = boost::mpl::deref< typename boost::mpl::max_element<type>::type  >::type::value;
    static const int size = boost::mpl::size<V>::value;
    static_assert( (ma - mi +1 == size), "Indices non contiguous in memory"); 
   }; 
  };

  // map this transform to the GI. GI2 is a boost::mpl::vector of boost::mpl::vector is memory position
  typedef typename boost::mpl::transform <MIndexList, one_index_group_to_memory_pos >::type  GI2;

  // get the first representant of each vector
  struct get_head { template<typename V> struct apply : boost::mpl::at_c<V,0> {}; };
  typedef typename boost::mpl::transform <GI2, get_head >::type  GI3;

  // now we need to count for each block how many other block have a (first) index smaller than their first index
  // !! count_if return an unsigned int, while the permutations metacode uses boost::mpl::int_ ---> breaks inversion !
  struct get_pos { template<typename V> struct apply : boost::mpl::int_<boost::mpl::count_if < GI3, boost::mpl::less< boost::mpl::_, V> >::value> {}; };
  typedef typename boost::mpl::transform <GI3, get_pos>::type new_memory_pos;

  template< typename V1> struct my_permutation : Permutations::perm_tag {typedef V1 V;};
  typedef my_permutation<new_memory_pos> P;
  typedef array_view<typename A::value_type, new_dim, Option::options< Option::memory_order_p<P> > > type;  // rest of options ? forward them ??

  struct _rt_data {
   mini_vector<size_t, type::rank> l; mini_vector<std::ptrdiff_t, type::rank> s;
   size_t i; A const & a; 
   _rt_data(A const & a_): i(0),a(a_) { for (size_t u=0; u<type::rank; ++u) {l[u]=1; s[u]=0;} }
  };

  struct _rt2 {
   _rt_data & d; int i;
   _rt2(_rt_data & d_): d(d_), i(0) {}
   template<typename I> void operator()(I) { 
    d.l[d.i] *= d.a.indexmap().lengths()[I::value]; 
    d.s[d.i] = (i==0 ? d.a.indexmap().strides()[I::value]: std::min(d.s[d.i], d.a.indexmap().strides()[I::value]));
    ++i;
   } 
  };

  struct _rt1 { 
   _rt_data & d;
   _rt1(_rt_data & d_): d(d_) {}
   template<typename V> void operator()(V) { _rt2 r(d); boost::mpl::for_each<typename V::type>(r); ++d.i; }
  };

  static type invoke(A const & a) {
   _rt_data r(a);
   boost::mpl::for_each<MIndexList>( _rt1(r) ); 
   typename type::indexmap_type im(r.l,r.s,0);
   return type(im,a.storage());
  }
 };

 template <int f, BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT( BOOST_PP_DEC(ARRAY_NRANK_MAX), int n, -1) > 
  struct m_index { 
   typedef boost::mpl::vector_c<int, f, BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(ARRAY_NRANK_MAX),n) > type; };  

#define AUX(z, N, unused) \
 template <BOOST_PP_ENUM_PARAMS(N,int n)> \
 struct m_index< BOOST_PP_ENUM_PARAMS(N,n) BOOST_PP_ENUM_TRAILING_PARAMS(BOOST_PP_SUB(ARRAY_NRANK_MAX, N), -1 BOOST_PP_INTERCEPT) > {\
  typedef boost::mpl::vector_c<int, BOOST_PP_ENUM_PARAMS(N,n) > type; };
 BOOST_PP_REPEAT_FROM_TO(1,ARRAY_NRANK_MAX, AUX, nil)
#undef AUX

#define AUX(z, N, unused) \
  template<BOOST_PP_ENUM_PARAMS(N,typename MI) , typename A >\
  typename group_indices_impl<A,boost::mpl::vector<BOOST_PP_ENUM_PARAMS(N,MI) > >::type\
  group_indices(A const & a,BOOST_PP_ENUM_PARAMS(N,MI) ) { return group_indices_impl<A,boost::mpl::vector< BOOST_PP_ENUM_PARAMS(N,MI) > >::invoke(a);}
  BOOST_PP_REPEAT_FROM_TO(1,ARRAY_NRANK_MAX, AUX, nil)
#undef AUX

}}//namespace triqs::arrays 
#endif
