
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

#ifndef TRIQS_ARRAYS_MATRIX_CACHE_H
#define TRIQS_ARRAYS_MATRIX_CACHE_H
#include "../matrix.hpp"
//#include "../utility/typeid_name.hpp"
#include <boost/scoped_ptr.hpp>
#include <boost/mpl/if.hpp>
/*
 * matrix_cache <A , Policy, ValueType, Opt>
 *
 * Given an object a of type A modeling HasMatrixImmutableInterface (e.g. matrix, matrix_view, matrix_expression) 
 * it provides a cache : upon calling, it returns a exposed_type (matrix &/const&, or matrix_view &/const &).
 * Policy can be  :
 *   - const_copy : the cache ensures that you have access to a const version of the data (avoids copy if possible)
 *   - copy : ensures you have a non-const copy : you can modify the copy BUT NOT the original data
 *   - reflexive : you have a copy that you can modify, and the data in back copied to a at the destruction of the cache object.
 * Cache is guaranteed to be contiguous.
 * Upon calling, it returns exposed_type either a matrix<ValueType, Opt> (const) & or a matrix_view<ValueType,Opt> (const) &
 * const in the case where Policy is const_copy
 * NB : No guarantee is made if A is modified between contruction and destruction of the cache.
 *      No guarantee is made on the validity of the ref to A, it is the caller responsability 
 *         (for optimisation, for a simple matrix, the cache just holds a & to the object).
 */
namespace triqs { namespace arrays { 

 namespace matrix_cache_policy { 
 enum policy_type { const_copy, copy, reflexive};
 }
 namespace matrix_cache_details { 
  using namespace matrix_cache_policy;
  enum internal_policy_type { always_copy, never_copy, may_copy};

  template<class A, policy_type P, class T, class Opt > struct _pol                    {static const internal_policy_type p = always_copy;};		
  template<class T, policy_type P, class Opt> struct _pol< matrix <T, Opt>,P,T,Opt>    {static const internal_policy_type p = never_copy;}; 
  template<class T, class Opt> struct _pol< matrix <T, Opt>,copy,T,Opt>          {static const internal_policy_type p = always_copy;}; 
  template<class T, policy_type P, class Opt> struct _pol< matrix_view<T,Opt>,P,T,Opt> {static const internal_policy_type p = may_copy;}; 
  template<class T, class Opt> struct _pol< matrix_view<T,Opt>,copy,T,Opt>       {static const internal_policy_type p = always_copy;}; 

  template <class A> struct compute_opt { typedef Option::Default type;};
  template <class T, class Opt>  struct compute_opt< matrix<T,Opt> > { typedef Opt type;};
  template <class T, class Opt>  struct compute_opt< matrix_view<T,Opt> > { typedef Opt type;};

  template< class A, policy_type Policy, class T = typename A::value_type, class Opt = typename compute_opt<A>::type, 
   internal_policy_type policy = _pol<A,Policy ,T,Opt>::p // impl only, don't use this 
    > class matrix_cache;

  template<class A, policy_type P, class T, class Opt> 
   class matrix_cache<A,P,T,Opt,never_copy> {
   typedef typename mpl::if_c<(P!=const_copy), A & , A const &>::type A_ref; 
   A_ref a;
   public : 
   typedef A_ref constructor_arg_type; 
   explicit matrix_cache( A_ref a_): a(a_) {}
   typedef A_ref exposed_type;
   exposed_type operator()() const { return a;}
   };

  template<policy_type P, class A, class B> struct back_copy { static void invoke(A a,B b) {} };
  template<class A, class B> struct back_copy<reflexive,A,B> { static void invoke (A & a ,B  & b) {triqs::arrays::make_view(a) = b;} };

  template< class A, policy_type P, class T, class Opt > 
   class matrix_cache<A,P,T,Opt,always_copy> {
    typedef typename mpl::if_c<(P==reflexive), A & , A const &>::type A_ref; 
    A_ref orig;
    typedef matrix<T,Opt> data_type;
    boost::scoped_ptr< data_type > data;
    bool init;
    void prep() { data.reset(new data_type(orig)); init = true; } 
    public : 
    typedef A_ref constructor_arg_type; 
    explicit matrix_cache( A_ref a_): orig(a_), init(false)  {}
    matrix_cache(matrix_cache const &c):orig(c.orig),init(false) {}
    ~matrix_cache() { if (init) back_copy<P,A_ref,data_type>::invoke(orig,*data);} // copy data back
    typedef typename mpl::if_c<(P!=const_copy),data_type &, data_type const &>::type exposed_type; 
    exposed_type operator()()  { if(!init) prep(); return *data;}
   };

  template< class ViewType, policy_type P, class T, class Opt > 
   class matrix_cache<ViewType,P,T,Opt,may_copy> {
    typedef matrix<T,Opt> data_type;
    boost::scoped_ptr< data_type  > data;
    boost::scoped_ptr< matrix_view<T,Opt>  > V_ptr;
    bool need_copy, init;
    ViewType orig; 
    matrix_view<T,Opt> * V;
    void prep(){data.reset(new data_type(*V)); V_ptr.reset(new matrix_view<T,Opt>(*data)); V=V_ptr.get(); init=true; } 
    public :
    typedef ViewType const & constructor_arg_type; 
    explicit matrix_cache (ViewType const & x): need_copy (!( x.indexmap().is_contiguous())), init(false), orig(x), V(&orig){
     if (need_copy) std::cerr<< " I need a copy "<< P <<std::endl;
    }
    matrix_cache(matrix_cache const &c):need_copy(c.need_copy),init(false),orig(c.orig),V(&orig) {}
    ~matrix_cache() { if (need_copy && init) back_copy<P,ViewType,data_type>::invoke(orig ,*(data));} // copy data back
    typedef typename mpl::if_c<(P!=const_copy),matrix_view<T,Opt>  &, matrix_view<T,Opt>  const &>::type exposed_type; 
    exposed_type operator()() { if (need_copy && (!init)) prep(); return *(V);}
   };

 }
 using matrix_cache_details::matrix_cache; 

}}//namespace triqs::arrays

#endif

