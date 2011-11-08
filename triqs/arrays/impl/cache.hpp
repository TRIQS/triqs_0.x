
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

#ifndef TRIQS_ARRAYS_ARRAY_CACHES_H
#define TRIQS_ARRAYS_ARRAY_CACHES_H
#include "../array.hpp"
#include "../matrix.hpp"
//#include "../tools/print_typeid.hpp"
/*
 * Purpose : 
 *  - Given an array/matrix,... or a view A, it provides a temporary copy, 
 *    reordered in C/Fortran order and contiguous.
 *  - If "reflexive", the data in back copied to A at destruction.
 *  - If A is already in the right order and compact, then just take a view.
 *
 *  - No guarantee is made if A is modified between contruction and destruction of the cache.
 *    array_cache is intended for local use only.
 *
 *
 * Usage : 
 *  - technical class, e.g. for HDF5 serialization.
 *  Just code the standard C ordered compact array, and 
 *  use the array_cache in the other cases.
 *  - When performance is not an issue.
 *  - cache < Reflexive , Order > (A)
 *  - the user functions are in the details:: namespace.
 *
 * Example : 
 *
 *  - details::cache_CC (A)
 * 
 */
namespace triqs { namespace arrays {
 namespace details_cache { 

  // cache_object_type : A,OrderTag ---> the type to use to make the cache. 
  template<typename A, typename OrderTag> 
   struct cache_object_type { typedef array<typename A::value_type,A::domain_type::rank,Option::options<OrderTag > > type;}; // default

  template<typename V, int rank, typename Opt, typename OrderTag>
   struct cache_object_type< array<V,rank,Opt>, OrderTag > { typedef  array<V,rank,Option::options<OrderTag> > type; };

  template<typename V, int rank, typename Opt, typename OrderTag>
   struct cache_object_type< array_view<V,rank,Opt>, OrderTag > { typedef  array<V,rank,Option::options<OrderTag> > type; };

  template <typename V, typename Opt, typename OrderTag>
   struct cache_object_type< matrix <V,Opt>, OrderTag > { typedef  matrix<V,Option::options<OrderTag> > type; };

  template <typename V, typename Opt, typename OrderTag>
   struct cache_object_type< matrix_view <V,Opt>, OrderTag > { typedef  matrix<V,Option::options<OrderTag> > type; };

  // -------------------------------------------------

  template<typename A, typename CacheObjectType, bool Reflexive> struct cache_impl;

  // the non reflexive case
  template<typename A, typename CacheObjectType> struct cache_impl<A, CacheObjectType, false> { 
   typedef typename CacheObjectType::view_type V;
   boost::shared_ptr<CacheObjectType> tmp;
   V v;
   cache_impl(A const & x): tmp(new CacheObjectType(x)), v(*tmp) {
#ifdef TRIQS_ARRAYS_CACHE_VERBOSE
    std::cerr<<"Cache : copy made" <<std::endl;
#endif
   }
   operator V const & () { return v;}
   V const & view() const { return v;}
  };

  // reflexive case : just add the back copy of the data at destruction
  template<typename A, typename CacheObjectType> 
   struct cache_impl<A, CacheObjectType, true> : cache_impl<A, CacheObjectType, false> { 
    typedef cache_impl<A, CacheObjectType, false> BaseType;
    A & a;
    cache_impl(A & x): BaseType(x),a(x) {}
    cache_impl(cache_impl const & X): BaseType(X),a(X.a){}
    ~cache_impl() { triqs::arrays::make_view(a) = this->v; } // copy data back
    operator typename BaseType::V &  () { return this->v;}
    typename BaseType::V & view() { return this->v;}
   };

  // special case when the objects are the same. No copy is made, except if data are not compact 
  template<typename A> struct cache_impl<A, A, false> { 
   typedef A CacheObjectType;
   typedef typename CacheObjectType::view_type V;
   A const & a;
   bool no_copy;
   boost::shared_ptr<CacheObjectType> tmp;
   V v;
   cache_impl(A const & x):a(x), 
   no_copy( x.indexmap().is_contiguous() ),
   tmp(no_copy ? NULL : new CacheObjectType(x)),
   v(no_copy ? x : *tmp) {
#ifdef TRIQS_ARRAYS_CACHE_VERBOSE
    std::cerr<<(no_copy ? "Cache : no copy made" : "Cache : type identical but copy made")<<std::endl;
#endif
   }
   operator V const & () { return v;}
   V const & view() const { return v;}
  };

  // reflexive case : just add the back copy of the data at destruction
  template<typename A> struct cache_impl<A, A, true> : cache_impl<A, A, false> { 
   typedef cache_impl<A, A, false> BaseType;
   A & a;
   cache_impl(A & x): BaseType(x),a(x) {}
   ~cache_impl() { if (!this->no_copy) triqs::arrays::make_view(a) = this->view(); } // copy data back
   operator typename BaseType::V &  () { return this->v;}
   typename BaseType::V & view() { return this->v;}
  };

  // the user routines
  template<typename A>
   cache_impl<A, typename cache_object_type<A,Tag::C>::type ,false>
   cache_CC(A const & a) { return cache_impl<A, typename cache_object_type<A,Tag::C>::type ,false>(a);}

  template<typename A>
   cache_impl<A, typename cache_object_type<A,Tag::C>::type ,true>
   cache_CR(A & a) { return cache_impl<A, typename cache_object_type<A,Tag::C>::type , true>(a);}

  template<typename A>
   cache_impl<A, typename cache_object_type<A,Tag::Fortran>::type ,false>
   cache_FC(A const & a) { return cache_impl<A, typename cache_object_type<A,Tag::Fortran>::type ,false>(a);}

  template<typename A>
   cache_impl<A, typename cache_object_type<A,Tag::Fortran>::type ,true>
   cache_FR(A & a) { return cache_impl<A, typename cache_object_type<A,Tag::Fortran>::type , true>(a);}

 }// details_cache

 namespace result_of { 
  template<bool Reflexive, typename OrderTag, typename A> 
   struct cache { typedef details_cache::cache_impl<A, typename details_cache::cache_object_type<A,OrderTag>::type , Reflexive> type; };
 }

 namespace details { 
  // put the user routine the right namespace
  using details_cache::cache_CC;
  using details_cache::cache_CR;
  using details_cache::cache_FC;
  using details_cache::cache_FR;
 }

}}//namespace triqs::arrays 
#endif

