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
#include <boost/shared_ptr.hpp>

namespace triqs { namespace arrays {

 template<typename A1, typename A2, typename Enable =void > struct need_copy_ct : mpl::true_{};

 template<typename A1, typename A2> struct need_copy_ct<A1,A2, typename boost::enable_if<is_value_or_view_class<A1> >::type> : 
  mpl::not_<indexmaps::IndexOrder::same_order<typename A1::indexmap_type::index_order_type, typename A2::indexmap_type::index_order_type> >{};

 template<typename DataType, typename CacheType, bool ct_need_copy = need_copy_ct<DataType,CacheType>::value > class const_cache;
 template<typename DataType, typename CacheType, bool ct_need_copy = need_copy_ct<DataType,CacheType>::value > class cache;

 template<typename D, typename C> const_cache<D,C>
  make_const_cache( D const & x) { return const_cache<D,C>(x);}

 template<typename D> const_cache<D, array<typename D::value_type, D::domain_type::rank,Option::C> >
  make_const_cache_C_order( D const & x) { return const_cache<D, array<typename D::value_type, D::domain_type::rank,Option::C> >(x);}

 template<typename D> cache<D, array<typename D::value_type, D::domain_type::rank,Option::C> >
  make_cache_C_order( D const & x) { return cache<D, array<typename D::value_type, D::domain_type::rank,Option::C> >(x);}

 template<typename D> const_cache<D, array<typename D::value_type, D::domain_type::rank,Option::Fortran> >
  make_const_cache_fortran_order( D const & x) { return const_cache<D, array<typename D::value_type, D::domain_type::rank,Option::Fortran> >(x);}

 template<typename D> cache<D, array<typename D::value_type, D::domain_type::rank,Option::Fortran> > 
  make_cache_fortran_order( D const & x) { return cache<D, array<typename D::value_type, D::domain_type::rank,Option::Fortran> >(x);}

 // ----------------- implementation  ----------------------------------

 template<typename A1, typename Enable =void > struct get_orig { typedef A1 type;};
 template<typename A1> struct get_orig<A1,typename boost::enable_if<is_value_or_view_class<A1> >::type> { typedef typename A1::view_type type;};

 // first case : copy is mandatory from compile time decision
 template<typename DataType, typename CacheType> 
  class const_cache<DataType,CacheType,true> {
   protected: 
   typedef typename get_orig<DataType>::type original_view_type;// view type if a regular class, regular type if an expression
   original_view_type  original_view;
   typedef typename CacheType::view_type final_view_type;
   mutable boost::shared_ptr< CacheType > copy;
   mutable boost::shared_ptr< final_view_type > view_ptr;
   mutable bool init;
   final_view_type * prep() const { 
    if (!init) {
#ifdef TRIQS_ARRAYS_CACHE_COPY_VERBOSE
     std::cerr<< " Cache : copy made "<<std::endl;
#endif
     copy.reset(new CacheType(original_view)); view_ptr.reset(new final_view_type(*copy)); init=true; 
    } 
    return view_ptr.get();
   }
   public :
   explicit const_cache (DataType const & x): init(false), original_view (x){}
   void update() { if (this->init) *(this->view_ptr) = this->original_view;} 
   operator final_view_type const & () const { return view();}
   final_view_type const & view () const { return *(prep());}
  };

 // second case : copy is NOT mandatory from compile time decision. Hence original_view can a priori be made a final_view_type
 template<typename DataType, typename CacheType> 
  class const_cache<DataType,CacheType,false> : public const_cache<DataType,CacheType,true> { 
   typedef typename CacheType::view_type final_view_type;
   public :
   const bool need_copy;
   explicit const_cache (DataType const & x): const_cache<DataType,CacheType,true> (x), need_copy (!(x.indexmap().is_contiguous())){ }
   void update()  { if (this->need_copy && this->init) *(this->view_ptr) = this->original_view;} 
   final_view_type const & view  () const {if (!need_copy) return this->original_view; else return *(this->prep());}
   operator final_view_type const & () const {return view();}
  };

 // Non const case : just add the back copy in the destructor 
 template<typename DataType, typename CacheType> 
  class cache<DataType,CacheType,true> : const_cache<DataType,CacheType,true> { 
   static_assert( is_value_or_view_class<DataType>::value, "non const cache only for regular classes and views, not expressions");
   typedef typename CacheType::view_type final_view_type;
   public :
   explicit cache (DataType const & x): const_cache<DataType,CacheType,true>  (x) {}
   ~cache() { back_update(); }
   void back_update() { if (this->init) this->original_view  = *(this->view_ptr);} 
   final_view_type &       view ()       {return *(this->prep());}
   final_view_type const & view () const {return *(this->prep());}
   operator final_view_type & ()             {return view();}
   operator final_view_type const & () const {return view();}
  };

 // Non const case : just add the back copy in the destructor 
 template<typename DataType, typename CacheType> 
  class cache<DataType,CacheType,false> : const_cache<DataType,CacheType,false> { 
   static_assert( is_value_or_view_class<DataType>::value, "non const cache only for regular classes and views, not expressions");
   typedef typename CacheType::view_type final_view_type;
   public :
   explicit cache (DataType const & x): const_cache<DataType,CacheType,false>  (x) {}
   ~cache() { back_update(); }
   void back_update() { if (this->need_copy && this->init) this->original_view  = *(this->view_ptr);} 
   final_view_type &       view ()       {if (! this->need_copy) return this->original_view; return *(this->prep());}
   final_view_type const & view () const {if (! this->need_copy) return this->original_view; return *(this->prep());}
   operator final_view_type & ()             {return view();}
   operator final_view_type const & () const {return view();}
  };

}}//namespace triqs::arrays 
#endif


