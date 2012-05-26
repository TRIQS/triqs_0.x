
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
#ifndef TRIQS_ARRAYS_COMPOUND_ASSIGN2_H_
#define TRIQS_ARRAYS_COMPOUND_ASSIGN2_H_
#include "./assignment.hpp"

#define TRIQS_DEFINE_COMPOUND_OPERATORS(MYTYPE)\
  template<typename RHS> MYTYPE & operator +=(RHS const & rhs) { compound_assignment<MYTYPE,RHS,'A'>(*this,rhs); return *this;  }\
  template<typename RHS> MYTYPE & operator -=(RHS const & rhs) { compound_assignment<MYTYPE,RHS,'S'>(*this,rhs); return *this;  }\
  template<typename RHS> MYTYPE & operator *=(RHS const & rhs) { compound_assignment<MYTYPE,RHS,'M'>(*this,rhs); return *this;  }\
  template<typename RHS> MYTYPE & operator /=(RHS const & rhs) { compound_assignment<MYTYPE,RHS,'D'>(*this,rhs); return *this;  }

namespace triqs { namespace arrays { 

 namespace details { template<typename LHS, typename RHS, char OP, typename Enable = void>  struct comp_assign_impl; }

 // puts the contents of RHS into LHS. LHS must be an indexmap_storage_pair
 // it is specialized in various cases for optimisation.
 template<typename LHS, typename RHS, char OP> 
  void compound_assignment  (LHS & lhs, const RHS & rhs)  { details::comp_assign_impl<LHS,RHS,OP>(lhs,rhs).invoke();}

 // -------- IMPLEMENTATION ----------------------------
 template <typename T> struct is_vector_or_view;

 namespace details { 

  template<typename A,typename B, char OP> struct _ops_;
  // faux : use make_expr and typededuction .....
  template<typename A,typename B> struct _ops_ <A,B,'A'> { static void invoke1 (A & a, B const & b) { a+=b;} };
  template<typename A,typename B> struct _ops_ <A,B,'S'> { static void invoke1 (A & a, B const & b) { a-=b;} };
  template<typename A,typename B> struct _ops_ <A,B,'M'> { static void invoke1 (A & a, B const & b) { a*=b;} };
  template<typename A,typename B> struct _ops_ <A,B,'D'> { static void invoke1 (A & a, B const & b) { a/=b;} };

  template<char OP> struct _is_AS: boost::mpl::false_{};
  template<> struct _is_AS<'A'>: boost::mpl::true_{};
  template<> struct _is_AS<'S'>: boost::mpl::true_{};

  template<char OP> struct _is_MD: boost::mpl::false_{};
  template<> struct _is_MD<'M'>: boost::mpl::true_{};
  template<> struct _is_MD<'D'>: boost::mpl::true_{};

  using indexmaps::foreach;

  // standard assignment for indexmap_storage_pair
  template<typename LHS, typename RHS, char OP> 
   struct comp_assign_impl<LHS,RHS, OP, typename boost::enable_if<boost::mpl::and_<is_isp<RHS,LHS>, boost::mpl::not_<is_vector_or_view<LHS> > > >::type> :
   _ops_<typename LHS::value_type, typename RHS::value_type,OP> { 
    typedef typename LHS::value_type value_type;
    typedef typename LHS::indexmap_type indexmap_type;
    typedef typename indexmap_type::domain_type::index_value_type index_value_type;
    LHS & lhs; const RHS & rhs;
    comp_assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}

    template<typename S1, typename S2>
     static typename boost::enable_if<storages::raw_copy_possible<S1,S2>, bool >::type
     rcp_impl(S1 & s1, S2 const & s2) {s1.raw_copy_from(s2); return true;}

    template<typename S1, typename S2>
     static typename boost::disable_if<storages::raw_copy_possible<S1,S2>, bool >::type
     rcp_impl(S1 & s1, S2 const & s2) { return false;} 

    void operator()(value_type & p, index_value_type const & key) { this->invoke1( p, rhs[key]);}
    void invoke ()  {
#ifdef TRIQS_ARRAYS_DEBUG
     if (!indexmaps::compatible_for_assignment(lhs.indexmap(), rhs.indexmap())) throw "Size mismatch";
#endif
     if ((indexmaps::raw_copy_possible(lhs.indexmap(), rhs.indexmap()))  &&
       (lhs.storage().size()== rhs.storage().size()) && rcp_impl(lhs.storage(),rhs.storage()) ) {} // SUPPRESS THIS  !!!!! 
     else {
#ifndef TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH 
      typename RHS::const_iterator it_rhs = rhs.begin();
      typedef typename RHS::const_iterator::indexmap_iterator_type RHS_mapit;
      typedef typename indexmaps::indexmap_iterator_adapter< RHS_mapit, typename LHS::indexmap_type >::type IT;
      iterator_adapter<false, IT, typename LHS::storage_type > it_lhs(lhs.indexmap(),lhs.storage());
      for (;it_lhs; ++it_lhs, ++it_rhs) { assert(it_rhs); invoke1(*it_lhs,*it_rhs); }
#else
      foreach(*this,lhs); 
#endif
     }
    }
   };

  // assignment for expressions RHS
  template<typename LHS, typename RHS, char OP> 
   struct comp_assign_impl<LHS,RHS,OP, typename boost::enable_if<boost::mpl::and_< ImmutableArray<RHS>, 
   boost::mpl::not_<Tag::check<Tag::has_special_infix<OP>,RHS> >,
   boost::mpl::not_< is_scalar_for<RHS,LHS > >,
   boost::mpl::not_< is_isp<RHS,LHS> > > >::type > : //is_expression<RHS> >::type > {    
    _ops_<typename LHS::value_type, typename RHS::value_type,OP> {    
    typedef typename LHS::value_type value_type;
    typedef typename LHS::indexmap_type indexmap_type;
    typedef typename indexmap_type::domain_type::index_value_type index_value_type;
    LHS & lhs; const RHS & rhs;
    comp_assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
    void operator()(value_type & p, index_value_type const & key) const { this->invoke1( p, rhs[key]);}
    void invoke() { 
#ifdef TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH 
     foreach(*this,lhs); 
#else
     typename LHS::storage_type & S(lhs.storage());
     for (typename LHS::indexmap_type::iterator it(lhs.indexmap());it; ++it)  invoke1( S[*it] ,rhs[it.indices()]) ;  
#endif
    }
   };

  // *= and /= for scalar RHS, except for vector (below).
  template<typename LHS, typename RHS, char OP> 
   struct comp_assign_impl<LHS,RHS,OP, 
   typename boost::enable_if<boost::mpl::and_<_is_MD<OP>, is_scalar_for<RHS,LHS >,  boost::mpl::not_<is_vector_or_view<LHS> > > >::type > :
    _ops_<typename LHS::value_type, RHS,OP> {    
     typedef typename LHS::value_type value_type;
     typedef typename LHS::indexmap_type::domain_type::index_value_type index_value_type;
     LHS & lhs; const RHS & rhs;
     comp_assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
     void operator()(value_type & p, index_value_type const & key) const { this->invoke1( p , rhs); }
     void invoke() { 
#ifdef TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH 
      foreach(*this,lhs);  // if contiguous : plain loop else foreach...
#else
      typename LHS::storage_type & S(lhs.storage());
      for (typename LHS::indexmap_type::iterator it(lhs.indexmap());it; ++it)  invoke1(S[*it] , rhs);  
#endif
     }
    };

  template<typename LHS, typename RHS, char OP> struct ca_special_;
  template<typename LHS, typename RHS> struct ca_special_< LHS,RHS,'A'> { static void invoke (LHS & lhs, RHS const & rhs) {rhs.assign_add_invoke(lhs);}};
  template<typename LHS, typename RHS> struct ca_special_< LHS,RHS,'S'> { static void invoke (LHS & lhs, RHS const & rhs) {rhs.assign_sub_invoke(lhs);}};
  template<typename LHS, typename RHS> struct ca_special_< LHS,RHS,'M'> { static void invoke (LHS & lhs, RHS const & rhs) {rhs.assign_mul_invoke(lhs);}};
  template<typename LHS, typename RHS> struct ca_special_< LHS,RHS,'D'> { static void invoke (LHS & lhs, RHS const & rhs) {rhs.assign_div_invoke(lhs);}};

  // assignment for anything that has a special (optimized) method for assignment
  template<typename LHS, typename RHS, char OP> 
   struct comp_assign_impl<LHS,RHS,OP,typename boost::enable_if<Tag::check<Tag::has_special_infix<OP>,RHS> >::type > { 
    LHS & lhs; const RHS & rhs;
    comp_assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
    void invoke ()  {  ca_special_< LHS,RHS,OP>::invoke(lhs,rhs); }
   };

 }// details
}}//namespace triqs::arrays
#endif

