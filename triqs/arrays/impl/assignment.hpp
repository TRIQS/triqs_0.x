
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

#ifndef TRIQS_ARRAYS_ASSIGN2_H_
#define TRIQS_ARRAYS_ASSIGN2_H_
#include "iterator_adapter.hpp"
#include "../indexmaps/cuboid/foreach.hpp"

// two ways of doing things... optimal one depends on compiler !
#define TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH

namespace triqs { namespace arrays { 

 namespace details { template<typename LHS, typename RHS, typename Enable = void>  struct assign_impl; }

 // puts the contents of RHS into LHS. LHS must be an indexmap_storage_pair
 // it is specialized in various cases for optimisation.
 template<typename LHS, typename RHS> 
  void assignment (LHS & lhs, const RHS & rhs)  { details::assign_impl<LHS,RHS>(lhs,rhs).invoke();}

 // -------- IMPLEMENTATION ----------------------------

//#define TRIQS_ARRAYS_TRACE_ASSIGN
 
 namespace details { 

  // RHS is considered to be an indexmap_storage_pair if it is one, ... except if it is the scalar type of hte LHS
  // think about an Array< Array<T,2> > e.g.
  template<class RHS,class LHS> struct is_isp : 
   boost::mpl::and_< Tag::check<Tag::indexmap_storage_pair,RHS>, boost::mpl::not_<is_scalar_for<RHS, LHS > > > {};

  // standard assignment for indexmap_storage_pair
  // for humans it means : 
  // if RHS and LHS are isp and NOT( LHS and RHS are both vectors)
  // in the both vector case, we have a specialized one in vector.hpp, using blas.
  template<typename LHS, typename RHS> 
   struct assign_impl<LHS,RHS, typename boost::enable_if<
   boost::mpl::and_< is_isp<RHS,LHS>, boost::mpl::not_<
   boost::mpl::and_<is_vector_or_view<LHS>,is_vector_or_view <RHS > > > > >::type> {
    typedef typename LHS::value_type value_type;
    typedef typename LHS::indexmap_type indexmap_type;
    typedef typename indexmap_type::domain_type::index_value_type index_value_type;
    LHS & lhs; const RHS & rhs;
    assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}

    template<typename S1, typename S2>
     static typename boost::enable_if<storages::raw_copy_possible<S1,S2>, bool >::type
     rcp_impl(S1 & s1, S2 const & s2) {s1.raw_copy_from(s2); return true;}

    template<typename S1, typename S2>
     static typename boost::disable_if<storages::raw_copy_possible<S1,S2>, bool >::type
     rcp_impl(S1 & s1, S2 const & s2) { return false;} 

    void operator()(value_type & p, index_value_type const & key) { 
     const_cast<typename boost::remove_const<value_type>::type &>(p)=rhs[key];}
    // hack : to be cleaned. When construction array<const T>(from expression, it does not work properly.
    // to fix more generally. Pb occurred in gf.

    void invoke ()  {
#ifdef TRIQS_ARRAYS_TRACE_ASSIGN
     std::cerr<<" assign for ISP"<<std::endl;
#endif
#ifdef TRIQS_ARRAYS_DEBUG
     if (!indexmaps::compatible_for_assignment(lhs.indexmap(), rhs.indexmap())) throw "Size mismatch";
#endif
     if ((indexmaps::raw_copy_possible(lhs.indexmap(), rhs.indexmap())) && 
       (lhs.storage().size()== rhs.storage().size()) && rcp_impl(lhs.storage(),rhs.storage()) ) {} 
     else {
#ifdef TRIQS_ARRAYS_TRACE_ASSIGN
     std::cerr<<" assign for ISP : no raw copy was possible"<<std::endl;
#endif
#ifndef TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH 
      typename RHS::const_iterator it_rhs = rhs.begin();
      typedef typename RHS::const_iterator::indexmap_iterator_type RHS_mapit;
      typedef typename indexmaps::indexmap_iterator_adapter< RHS_mapit, typename LHS::indexmap_type >::type IT;
      iterator_adapter<false, IT, typename LHS::storage_type > it_lhs(lhs.indexmap(),lhs.storage());
      for (;it_lhs; ++it_lhs, ++it_rhs) { assert(it_rhs); *it_lhs = *it_rhs; }
#else
      //std::cerr<<" vector = with foreach "<<std::endl;  
      indexmaps::foreach(*this,lhs); 
#endif
     }
    }
   };

  // assignment for expressions RHS
  template<typename LHS, typename RHS> 
   struct assign_impl<LHS,RHS,typename boost::enable_if<is_expression<RHS> >::type > { 
    typedef typename LHS::value_type value_type;
    typedef typename LHS::indexmap_type indexmap_type;
    typedef typename indexmap_type::domain_type::index_value_type index_value_type;
    LHS & lhs; const RHS & rhs;
    assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
    void operator()(value_type & p, index_value_type const & key) const { p=rhs[key];}
    void invoke() { 
#ifdef TRIQS_ARRAYS_TRACE_ASSIGN
     std::cerr<<" assign for expression"<<std::endl;
#endif
#ifdef TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH 
     indexmaps::foreach(*this,lhs); 
#else
     typename LHS::storage_type & S(lhs.storage());
     for (typename LHS::indexmap_type::iterator it(lhs.indexmap());it; ++it)  S[*it] = rhs[it.indices()];  
#endif
    } 
   };

  // assignment for scalar RHS // write a specific one if it is not a view : plain loop
  // beware : for matrix, assign to a scalar will make the matrix scalar, as it should
  // it will be specialized in the matrix class.... Cf matrix class
  template<typename LHS, typename RHS> 
   struct assign_impl<LHS,RHS,typename boost::enable_if<is_scalar_for<RHS,LHS > >::type > { 
    typedef typename LHS::value_type value_type;
    typedef typename LHS::indexmap_type indexmap_type;
    typedef typename indexmap_type::domain_type::index_value_type index_value_type;
    LHS & lhs; const RHS & rhs;
    assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
    void operator()(value_type & p, index_value_type const & key) const { p = rhs;}
    void invoke() {  
#ifdef TRIQS_ARRAYS_TRACE_ASSIGN
     std::cerr<<" scalar assign"<<std::endl;
#endif
#ifdef TRIQS_ARRAYS_ASSIGN_ISP_WITH_FOREACH 
     indexmaps::foreach(*this,lhs);  // if contiguous : plain loop else foreach...
#else
     typename LHS::storage_type & S(lhs.storage());
     for (typename LHS::indexmap_type::iterator it(lhs.indexmap());it; ++it)  S[*it] = rhs;  
#endif
    }
   };

  // assignment for anything that has a special (optimized) method for assignment
  template<typename LHS, typename RHS> 
   struct assign_impl<LHS,RHS,typename boost::enable_if< has_special_assign<RHS> >::type > { 
    LHS & lhs; const RHS & rhs;
    assign_impl(LHS & lhs_, const RHS & rhs_): lhs(lhs_), rhs(rhs_) {}
    void invoke ()  { 
#ifdef TRIQS_ARRAYS_TRACE_ASSIGN
     std::cerr<<" special assign invoke "<<std::endl;
#endif
     rhs.assign_invoke(lhs);
    }
   };

 }// details
}}//namespace triqs::arrays 
#endif

