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
#ifndef TRIQS_ARRAYS_SLIDING_VIEW_H
#define TRIQS_ARRAYS_SLIDING_VIEW_H
#include "./array.hpp"
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_trailing.hpp>
//#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

namespace triqs { namespace arrays {

 template<typename ArrayType,int Pos > class sliding_view;

#ifndef DO_NOT_DEFINE_ME
 // human version of the class, the preprocessor generalisation is next..
 template<typename ArrayType > class sliding_view<ArrayType,2> : TRIQS_MODEL_CONCEPT(MutableMatrix)  {
  typename ArrayType::view_type A; size_t n;
  public :
  typedef typename ArrayType::value_type value_type;
  sliding_view (ArrayType const & A_) : A(A_), n(0){}

  typedef indexmaps::slicer<typename ArrayType::indexmap_type , range , range,size_t,ellipsis> slicer_t;
  typedef typename slicer_t::r_type indexmap_type;
  indexmap_type indexmap() const { return slicer_t::invoke(A.indexmap() , range() , range(),n, ellipsis()); }
  typedef typename indexmap_type::domain_type domain_type;
  domain_type domain() const { return indexmap().domain();}
  typename ArrayType::storage_type const & storage() const { return A.storage();}
  typename ArrayType::storage_type & storage() { return A.storage();}
  size_t size() const { return A.len(2);}
  void set(size_t p) { n=p;}

  value_type const * restrict data_start() const { return A.data_start(); }// &storage_[indexmap_.start_shift()];}
  value_type * restrict data_start() { return A.data_start(); }//	eturn &storage_[indexmap_.start_shift()];}

  template<typename RHS> sliding_view & operator=(const RHS & X) {triqs_arrays_assign_delegation(*this,X); return *this; }
  TRIQS_DEFINE_COMPOUND_OPERATORS(sliding_view);

  template< typename A0 , typename A1 , typename ... Args> value_type & operator() ( A0 &&a0 , A1 &&a1 , Args && ... args)
  { return A( std::forward<A0>(a0) , std::forward<A1>(a1) , n,std::forward<Args>(args)...);}
  template< typename A0 , typename A1 , typename ... Args> value_type const & operator() ( A0 &&a0 , A1 &&a1 , Args && ... args) const
  { return A( std::forward<A0>(a0) , std::forward<A1>(a1) , n,std::forward<Args>(args)...);}
};
#else
#define AUX0(z,P,NNN) std::forward<A##P>(a##P),
#define AUX1(z,P,NNN) A##P && a##P,
#define TEXT(z, n, text) text
#define IMPL(z, POS, unused)\
 template<typename ArrayType >\
class sliding_view<ArrayType,POS> : TRIQS_MODEL_CONCEPT(ImmutableCuboidArray) {\
 typedef typename ArrayType::value_type value_type;\
 typename ArrayType::view_type A; size_t n;\
 public :\
	 sliding_view (ArrayType const & A_) : A(A_), n(0){}\
 \
 typedef indexmaps::slicer<typename ArrayType::indexmap_type BOOST_PP_ENUM_TRAILING(POS, TEXT, range),size_t,ellipsis> slicer_t;\
 typedef typename slicer_t::r_type indexmap_type;\
 indexmap_type indexmap() const { return slicer_t::invoke(A.indexmap() BOOST_PP_ENUM_TRAILING(POS, TEXT, range()),n, ellipsis()); }\
 typename ArrayType::storage_type storage() const { return A.storage();}\
 \
 size_t size() const { return A.len(POS);}\
 void set(size_t p) { n=p;}\
 \
 template<BOOST_PP_ENUM_PARAMS(POS,typename A) BOOST_PP_COMMA_IF(POS) typename ... Args>\
 value_type & operator() (BOOST_PP_REPEAT(POS,AUX1,nil) Args && ... args)\
 { return A(BOOST_PP_REPEAT(POS,AUX0,nil) n,std::forward<Args>(args)...);}\
 \
 template<BOOST_PP_ENUM_PARAMS(POS,typename A) BOOST_PP_COMMA_IF(POS) typename ... Args>\
 value_type const & operator() (BOOST_PP_REPEAT(POS,AUX1,nil)  Args && ... args) const \
 { return A(BOOST_PP_REPEAT(POS,AUX0,nil) n,std::forward<Args>(args)...);}\
};

BOOST_PP_REPEAT(ARRAY_NRANK_MAX , IMPL, nil);
#undef IMPL
#undef AUX0
#undef AUX1
#undef TEXT
#endif

template<int Pos, typename ArrayType>
sliding_view<ArrayType,Pos> make_sliding_view(ArrayType const & A) { return sliding_view<ArrayType,Pos> (A);}

/*
 * template<typename ArrayType>
 class sliding_view1 {
 typedef typename ArrayType::view_type array_view_t;
// what about Opt and To here ?
typedef array_view<typename ArrayType::value_type, ArrayType::rank-1>  exposed_array_view_t;
exposed_array_view_t V;
typename exposed_array_view_t::indexmap_type * vim;
//size_t l_max; // for checks
std::ptrdiff_t shift_, shift_0;
std::ptrdiff_t & start_shift_ref;
public :

template<typename ArrayType2, typename ... Args>
sliding_view1 (ArrayType2 const & A, size_t n, Args && ... args) :
V(A(std::forward<Args>(args)...)),
vim (const_cast<typename exposed_array_view_t::indexmap_type *>(&V.indexmap())),
start_shift_ref(const_cast<std::ptrdiff_t&>(V.indexmap().start_shift())) {
shift_ = A.indexmap().strides()[n];
shift_0 = vim->start_shift();

//l_max = A.indexmap().lengths()[n]
}

void reset() { vim->set_start_shift(shift_0);}

void shift (std::ptrdiff_t delta) { vim->set_start_shift(vim->start_shift() + delta * shift_);}

exposed_array_view_t & operator() (size_t n) { start_shift_ref = shift_0 + n * shift_; return V;}

void operator++() { start_shift_ref += shift_; } //vim->set_start_shift(vim->start_shift() + shift_);}


void operator--() { vim->set_start_shift(vim->start_shift() - shift_);}

exposed_array_view_t const & operator()() const { return V;}
exposed_array_view_t & operator()() { return V;}
};

template<typename ArrayType, typename ... Args>
sliding_view1<ArrayType> make_sliding_view1(ArrayType const & A,size_t n,Args && ... args) { return sliding_view1<ArrayType> (A,n,std::forward<Args>(args)...);}
*/

}}
#endif

