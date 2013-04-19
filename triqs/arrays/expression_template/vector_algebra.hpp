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
#ifndef TRIQS_ARRAYS_EXPRESSION_VECTOR_ALGEBRA_H
#define TRIQS_ARRAYS_EXPRESSION_VECTOR_ALGEBRA_H
#include "./tools.hpp"
#include "../matrix.hpp"
namespace triqs { namespace arrays {

 template<typename Tag, typename L, typename R> 
  struct vector_expr : TRIQS_MODEL_CONCEPT(ImmutableVector) { 
   typedef typename keeper_type<L>::type L_t;
   typedef typename keeper_type<R>::type R_t;
   L_t l; R_t r;
   template<typename LL, typename RR> vector_expr(LL && l_, RR && r_) : l(std::forward<LL>(l_)), r(std::forward<RR>(r_)) {}
   static_assert( get_rank<R_t>::value==0 || get_rank<L_t>::value==0 || get_rank<L_t>::value == get_rank<R_t>::value, "rank mismatch in array operations");
   typedef typename std::result_of<operation<Tag>(typename L_t::value_type,typename R_t::value_type)>::type  value_type;
   typedef typename std::remove_reference<typename std::result_of<combine_domain(L_t,R_t)>::type>::type domain_type;
 
   domain_type domain() const  { return combine_domain()(l,r); } 
   mini_vector<size_t,1> shape() const { return this->domain().lengths();} 
   //size_t dim0() const { return this->domain().lengths()[0];} 
   size_t size() const { return this->domain().lengths()[0];} 

   //template<typename KeyType> value_type operator[](KeyType && key) const { return operation<Tag>()(l[std::forward<KeyType>(key)] , r[std::forward<KeyType>(key)]);}
   template<typename ... Args> value_type operator()(Args && ... args) const { return operation<Tag>()(l(std::forward<Args>(args)...) , r(std::forward<Args>(args)...));}
   friend std::ostream &operator <<(std::ostream &sout, vector_expr const &expr){return sout << "("<<expr.l << " "<<operation<Tag>::name << " "<<expr.r<<")" ; }
  };

 template<typename L>  // a special case : the unary operator !
  struct vector_unary_m_expr : TRIQS_MODEL_CONCEPT(ImmutableVector) { 
   typedef typename keeper_type<L>::type L_t;
   L_t l; 
   template<typename LL> vector_unary_m_expr(LL && l_) : l(std::forward<LL>(l_)) {}
   typedef typename L_t::value_type value_type;
   typedef typename L_t::domain_type domain_type;

   domain_type domain() const  { return l.domain(); } 
   mini_vector<size_t,1> shape() const { return this->domain().lengths();} 
   //size_t dim0() const { return this->domain().lengths()[0];} 
   size_t size() const { return this->domain().lengths()[0];} 

   //template<typename KeyType> value_type operator[](KeyType&& key) const {return -l[key];} 
   template<typename ... Args> value_type operator()(Args && ... args) const { return -l(std::forward<Args>(args)...);}
   friend std::ostream &operator <<(std::ostream &sout, vector_unary_m_expr const &expr){return sout << '-'<<expr.l; }
  };

 // Now we can define all the C++ operators ...
#define DEFINE_OPERATOR(TAG, OP, TRAIT1, TRAIT2) \
 template<typename A1, typename A2>\
 typename std::enable_if<TRAIT1<A1>::value && TRAIT2 <A2>::value, vector_expr<tags::TAG, A1,A2>>::type\
 operator OP (A1 const & a1, A2 const & a2) { return vector_expr<tags::TAG, A1,A2>(a1,a2);} 
 DEFINE_OPERATOR(plus,       +, ImmutableVector,ImmutableVector);
 DEFINE_OPERATOR(minus,      -, ImmutableVector,ImmutableVector);
 DEFINE_OPERATOR(multiplies, *, is_in_ZRC,ImmutableVector);
 DEFINE_OPERATOR(multiplies, *, ImmutableVector,is_in_ZRC);
 DEFINE_OPERATOR(divides,    /, ImmutableVector,is_in_ZRC);
#undef DEFINE_OPERATOR

 // the unary is special
 template<typename A1> typename std::enable_if<ImmutableVector<A1>::value, vector_unary_m_expr<A1>>::type
  operator - (A1 const & a1) { return {a1};} 

 template<typename Expr > vector<typename Expr::value_type>
  make_vector( Expr const & e) { return vector<typename Expr::value_type>(e);}

}}//namespace triqs::arrays
#endif
