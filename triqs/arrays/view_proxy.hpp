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
#ifndef TRIQS_ARRAYS_VIEW_PROXY_H
#define TRIQS_ARRAYS_VIEW_PROXY_H
#include "./array.hpp"
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_trailing.hpp>
//#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

namespace triqs { namespace arrays {

 template<typename ArrayType,int Pos > class view_proxy;
 template<typename ArrayType,int Pos > class const_view_proxy;

 // to do : separate the array and the matrix case.
 // generalize with preprocessor (draft below)

 // write concept mutable down and clean it (dim0, dim1, len(i), ...) 
#ifndef DO_NOT_DEFINE_ME
 // human version of the class, the preprocessor generalisation is next..
 template<typename ArrayType > class const_view_proxy<ArrayType,2> : TRIQS_MODEL_CONCEPT(ImmutableMatrix)   {
  protected:
   ArrayType const * A; size_t n;
  public :
   typedef typename ArrayType::value_type value_type;
   const_view_proxy (ArrayType const &  A_, size_t n_=0) : A(&A_), n(n_){}
   typedef indexmaps::slicer<typename ArrayType::indexmap_type , range , range,size_t,ellipsis> slicer_t;
   typedef typename slicer_t::r_type indexmap_type;
   typedef typename indexmap_type::domain_type domain_type;
   indexmap_type indexmap() const { return slicer_t::invoke(A->indexmap() , range() , range(),n, ellipsis()); }
   domain_type domain() const { return indexmap().domain();}
   size_t len(int i) const { return A->len(i);}
   size_t dim0() const { return A->len(0);}
   size_t dim1() const { return A->len(1);}
   // shape ?
   typename ArrayType::storage_type const & storage() const { return A->storage();}
   template< typename A0 , typename A1 , typename ... Args> value_type const & operator() ( A0 &&a0 , A1 &&a1 , Args && ... args) const
   { return (*A)( std::forward<A0>(a0) , std::forward<A1>(a1) , n,std::forward<Args>(args)...);}
 };

 template<typename ArrayType > class view_proxy<ArrayType,2> : public const_view_proxy<ArrayType,2>, TRIQS_MODEL_CONCEPT(MutableMatrix){
  public :
  typedef const_view_proxy<ArrayType,2> B;
  view_proxy (ArrayType & A_, size_t n_=0) : B(A_,n_){}
  template<typename RHS> view_proxy & operator=(const RHS & X) {triqs_arrays_assign_delegation(*this,X); return *this; }
  TRIQS_DEFINE_COMPOUND_OPERATORS(view_proxy);
  template< typename A0 , typename A1 , typename ... Args> typename B::value_type & operator() ( A0 &&a0 , A1 &&a1 , Args && ... args) const 
   { return (*(const_cast<ArrayType*>(this->A)))( std::forward<A0>(a0) , std::forward<A1>(a1) , this->n,std::forward<Args>(args)...);}
 };

#ifdef RIQS_COMPLETE__
 template<typename ArrayType > class view_proxy<ArrayType,2> : TRIQS_MODEL_CONCEPT(MutableMatrix),TRIQS_MODEL_CONCEPT(MutableCuboidArray)   {
   ArrayType * A; size_t n;
  public :
  typedef typename ArrayType::value_type value_type;
  view_proxy (ArrayType & A_, size_t n_=0) : A(&A_), n(n_){}

  typedef indexmaps::slicer<typename ArrayType::indexmap_type , range , range,size_t,ellipsis> slicer_t;
  typedef typename slicer_t::r_type indexmap_type;
  typedef typename indexmap_type::domain_type domain_type;
  indexmap_type indexmap() const { return slicer_t::invoke(A->indexmap() , range() , range(),n, ellipsis()); }
  domain_type domain() const { return indexmap().domain();}
  size_t len(int i) const { return A->len(i);}
  size_t dim0() const { return A->len(0);}
  size_t dim1() const { return A->len(1);}
  
  typename ArrayType::storage_type const & storage() const { return A->storage();}

  template<typename RHS> view_proxy & operator=(const RHS & X) {triqs_arrays_assign_delegation(*this,X); return *this; }
  TRIQS_DEFINE_COMPOUND_OPERATORS(view_proxy);

  template< typename A0 , typename A1 , typename ... Args> value_type & operator() ( A0 &&a0 , A1 &&a1 , Args && ... args) const 
   { return (*A)( std::forward<A0>(a0) , std::forward<A1>(a1) , n,std::forward<Args>(args)...);}
 };
#endif

 #else
#define AUX0(z,P,NNN) std::forward<A##P>(a##P),
#define AUX1(z,P,NNN) A##P && a##P,
#define TEXT(z, n, text) text
#define IMPL(z, POS, unused)\
 template<typename ArrayType >\
 class view_proxy<ArrayType,POS> : TRIQS_MODEL_CONCEPT(ImmutableCuboidArray) {\
  typedef typename ArrayType::value_type value_type;\
  typename ArrayType::view_type A; size_t n;\
  public :\
	  view_proxy (ArrayType const & A_) : A(A_), n(0){}\
  \
  typedef indexmaps::slicer<typename ArrayType::indexmap_type BOOST_PP_ENUM_TRAILING(POS, TEXT, range),size_t,ellipsis> slicer_t;\
  typedef typename slicer_t::r_type indexmap_type;\
  indexmap_type indexmap() const { return slicer_t::invoke(A->indexmap() BOOST_PP_ENUM_TRAILING(POS, TEXT, range()),n, ellipsis()); }\
  typename ArrayType::storage_type storage() const { return A->storage();}\
  \
  size_t size() const { return A->len(POS);}\
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
  view_proxy<ArrayType,Pos> make_view_proxy(ArrayType const & A) { return view_proxy<ArrayType,Pos> (A);}

}}
#endif

