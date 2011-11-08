
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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

#ifndef TRIQS_ARRAYS_PERMUTATIONS_H
#define TRIQS_ARRAYS_PERMUTATIONS_H

#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/count.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/max_element.hpp>
#include <boost/mpl/min_element.hpp>
#include <boost/mpl/find.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp> 

#include <ostream>
#include <sstream>
#include <vector>

namespace triqs { namespace arrays { namespace Permutations { 
 struct perm_tag{};
 namespace mpl = boost::mpl;

#define NMAX ARRAY_NRANK_MAX

 /**
  *  permutation<n1,n2,n3> is the permutation [n1,n2,n3]
  *  e.g. permutation<2,3,1> is a 3 cycle
  */
#define TEXT(z, nval, text) int n##nval=-1
 template <int N , BOOST_PP_ENUM(BOOST_PP_DEC(NMAX), TEXT, dd)> 
#undef TEXT
  struct permutation: perm_tag { 
   typedef typename mpl::push_front< typename permutation< BOOST_PP_ENUM_PARAMS(BOOST_PP_DEC(NMAX), n) >::V , mpl::int_<N> >::type V;
   static_assert( (mpl::count<V, mpl::int_<N> >::type::value==1) ,  "Permutation : error, index appears twice");
  };
 template<> struct permutation<-1>: perm_tag { typedef mpl::vector<> V;};

 typedef permutation<-1> empty;

 template <int Head, typename P>
  struct push_front : perm_tag {typedef typename mpl::push_front< typename P::V, mpl::int_<Head> >::type V;}; 

 template <int Head, typename P>
  struct push_back  : perm_tag {typedef typename mpl::push_back < typename P::V, mpl::int_<Head> >::type V;}; 

 //  identity<N> is the identity permutation with N elements. [0, ..., N-1]
 template<int N> struct identity:push_back<N-1,identity<N-1> > {};
 template<> struct identity<0>:perm_tag {typedef mpl::vector<> V;};

 // reverse_identity<N> is [N-1, N-2, ... ,0] 
 template<int N> struct reverse_identity:push_front < N-1, reverse_identity<N-1> > {}; 
 template<> struct reverse_identity<0> : perm_tag  { typedef  mpl::vector<> V;};

 // is_permutation<P>::value == true iif P is a permutation. Unicity is by construction.
 template <typename P, int shift =0>
  struct is_permutation : mpl::bool_< ((mpl::max_element<P>::type::value ==  mpl::size<P>::type::value -1) &&  (mpl::min_element<P>::type::value ==  0) )> {};

 // length of the permutation
 template <typename P> struct length : mpl::size<typename P::V> {};

 // 
 template <typename P1, typename P2 >
  struct is_equal : mpl::equal<typename P1::V, typename P2::V> {};

 //  Transform into a vector
 struct to_vec_ { 
  std::vector<int> & V; to_vec_ (std::vector<int> & V_):V(V_){}
  template<typename T> void operator()(T) const { V.push_back(T::value);}
 };
 template <typename P> inline void to_vector (std::vector<int> & v) { mpl::for_each<typename P::V>(to_vec_(v));}

 //  Transform into a string
 template <typename P>
  inline std::string to_string () { 
   std::vector<int> V; to_vector<P>(V);
   std::stringstream s; s<<"[";
   for (unsigned int u=0; u<V.size(); ++u) s << (u>0 ? ", ": "") <<V[u]; 
   s<<"]"; return s.str();
  }

 // eval<P,N>::value is P(N). If Index is out of range, it generates a static_assert 
 template<typename P, int index>
  struct eval  { //: mpl::at<typename P::V,mpl::int_<index> > {};
   static_assert ( (index < mpl::size<typename P::V>::value), "Index out of range");
   typedef typename mpl::at<typename P::V,mpl::int_<index> >::type type;
   enum { value = type::value};
  };  

 /** 
  * permutations::find<y,P>::value is :
  *    -   X such that P(X)==y 
  *    -  (-1) if target is too large (or not found)    ?????
  *    */
 template < int y, typename P> struct find {
  static_assert ( ((y>=0) && (y < length<P>::value)), "Index out of range");
  typedef typename mpl::int_< mpl::find<typename P::V,mpl::int_<y> >::type::pos::value> type;
  enum { value = type::value};
 };

 // compose<P1,P2> is P1 * P2
 template <typename P1, typename P2, int c> struct comp_impl : push_back< eval<P1, eval<P2,c-1>::value>::value, comp_impl<P1,P2, c-1> > {};
 template <typename P1, typename P2> struct comp_impl<P1,P2,0>  { typedef  mpl::vector<> V;};
 template <typename P1, typename P2> struct compose : comp_impl<P1, P2, length<P1>::value >  {};

 // inverse<P> is the inverse permutation
 template<typename P, int c > struct inv_impl : push_back< find <c-1, P>::value, inv_impl<P, c-1> > {};
 template<typename P> struct inv_impl<P,0>  { typedef  mpl::vector<> V;};
 template<typename P> struct inverse : inv_impl<P, length<P>::value > {};

 /// Print a permutation 
 template <typename T> typename boost::enable_if<boost::is_base_of<perm_tag,T>, std::ostream &>::type
  operator<<(std::ostream & out, const T & p) { return out<<to_string<T>(); }

}

 using Permutations::permutation;

}}//namespace triqs::arrays 
#undef NMAX
#endif

