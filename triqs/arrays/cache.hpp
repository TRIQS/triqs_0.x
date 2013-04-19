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
#ifndef TRIQS_ARRAYS_CACHE_H
#define TRIQS_ARRAYS_CACHE_H
#include "./array.hpp"
#include <triqs/utility/view_tools.hpp>
namespace triqs { namespace arrays {

 template<typename DataType, typename CacheType> class const_cache;
 template<typename DataType, typename CacheType> class cache;
 
 template<typename D> 
  const_cache<D, array<typename D::value_type, D::domain_type::rank> > 
  make_const_cache( D const & x, memory_layout< D::domain_type::rank> ml = memory_layout< D::domain_type::rank>('C')  ) { 
   return const_cache<D, array<typename D::value_type, D::domain_type::rank> > (x,ml);
  }

 template<typename D> 
  cache<D, array<typename D::value_type, D::domain_type::rank> > 
  make_cache( D const & x, memory_layout< D::domain_type::rank> ml = memory_layout< D::domain_type::rank>('C') ) { 
   return cache<D, array<typename D::value_type, D::domain_type::rank> > (x,ml);
  }

 // ----------------- implementation  ----------------------------------

 // The type of A1, and A2 can already imply a copy. Compile time decision.
 template<typename A1, typename A2, typename Enable =void > struct need_copy_ct : mpl::true_{};
 template<typename A1, typename A2> struct need_copy_ct<A1,A2, ENABLE_IF(is_amv_value_or_view_class<A1>)> :  
  mpl::not_<std::is_same< typename A1::value_type, typename A2::value_type> > {};
  //mpl::not_<indexmaps::IndexOrder::same_order<typename A1::indexmap_type::index_order_type, typename A2::indexmap_type::index_order_type> >{};
 
 template<typename DataType, typename CacheType> class const_cache {
  protected:
   typedef memory_layout< DataType::domain_type::rank> ml_t;
   const ml_t ml;
   typename view_type_if_exists_else_type<DataType>::type keeper;
   static const bool CT_need_copy = need_copy_ct<DataType,CacheType>::value;
   const bool need_copy;
   typedef typename CacheType::view_type exposed_view_type;
   struct internal_data {
    CacheType copy;
    exposed_view_type view;
    internal_data(const_cache const & P,ml_t ml) : copy(CacheType(P.keeper,ml)), view(copy) {
#ifdef TRIQS_ARRAYS_CACHE_COPY_VERBOSE
     std::cerr<< " Cache : copy made "<< std::endl<< " -- TRACE = --" << std::endl << triqs::utility::stack_trace() << std::endl;
#endif
    }
   };
   friend struct internal_data;   
   mutable boost::shared_ptr<internal_data> _id;   
   internal_data  & id() const { if (!_id) _id= boost::make_shared<internal_data>(*this,ml); return *_id;}

   //exposed_view_type       & view2 ()       { if (need_copy) return id().view; else return keeper;}    
   //exposed_view_type const & view2 () const { if (need_copy) return id().view; else return keeper;} 

   // avoid compiling the transformation keeper-> exposed_view_type when it does not make sense
   exposed_view_type       & view1 (mpl::false_)       { if (need_copy) return id().view; else return keeper; }    
   exposed_view_type const & view1 (mpl::false_) const { if (need_copy) return id().view; else return keeper;} 

   exposed_view_type       & view1 (mpl::true_)       { return id().view; }     
   exposed_view_type const & view1 (mpl::true_) const { return id().view; } 

   exposed_view_type       & view2 ()       { return view1(mpl::bool_<CT_need_copy>());}
   exposed_view_type const & view2 () const { return view1(mpl::bool_<CT_need_copy>());}

   template<typename D> 
    typename std::enable_if< ! is_amv_value_or_view_class<D>::value, bool>::type 
    need_copy_dynamic ( D const & x) const { return false;}
   
   template<typename D> 
    typename std::enable_if<is_amv_value_or_view_class<D>::value, bool>::type 
    need_copy_dynamic ( D const & x) const { return (x.indexmap().memory_indices_layout() != ml );}

  public :
   const_cache (DataType const & x, ml_t ml_ = ml_t()): 
    ml(ml_),keeper (x), 
    need_copy ( CT_need_copy || need_copy_dynamic(x)  || (!has_contiguous_data(x)) ) {}
   void update() { if (need_copy && _id) id().view = keeper;} 
   exposed_view_type     const & view () const { return view2();}
   operator exposed_view_type const & () const { return view2();}
 };

 // Non const case : just add the back copy in the destructor 
 template<typename DataType, typename CacheType> class cache : const_cache<DataType,CacheType> { 
  static_assert( is_amv_value_or_view_class<DataType>::value, "non const cache only for regular classes and views, not expressions");
  typedef const_cache<DataType,CacheType> B;
  typedef typename B::ml_t ml_t;
  public :
  cache (DataType const & x, ml_t ml = ml_t() ): B(x,ml) {}
  ~cache() { back_update(); }
  void back_update() { if (this->need_copy) this->keeper = this->id().view;}
  typename B::exposed_view_type &       view ()           {return B::view2();}
  typename B::exposed_view_type const & view () const     {return B::view2();}
  operator typename B::exposed_view_type & ()             {return B::view2();}
  operator typename B::exposed_view_type const & () const {return B::view2();}
 };
}}//namespace triqs::arrays 
#endif


