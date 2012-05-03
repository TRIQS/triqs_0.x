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

#ifndef TRIQS_ARRAYS_QCACHE_H
#define TRIQS_ARRAYS_QCACHE_H
#include "./matrix.hpp"
#include "./vector.hpp"
#include <boost/scoped_ptr.hpp>
#include <boost/mpl/if.hpp>
namespace triqs { namespace arrays { 

 /*
  * const_qcache.
  * Given A, a matrix (or vector/array) it presents via the () operator
  *  - a const & to the matrix if A is a value class (matrix, vector, array, ..)
  *  - a const & to a new value class is A is a an expression 
  *  - a const & to a view if A is a view. If the view was not contiguous, it is a view to a 
  *    temporary regrouped value class.
  *
  *  const_qcache is NOT copyable. It should be used in local use only.
  *  [Rational : in simple case, like A =matrix, it allows to keep only a const &, which is quicker
  *
  */
 template<typename A, typename Enable=void> struct const_qcache;

 template<typename A> class const_qcache< A, typename boost::enable_if<is_value_class<A> >::type> : boost::noncopyable { 
  A const & x;
  public:
  const_qcache(A const & x_):x(x_){}
  typedef A const & exposed_type; 
  exposed_type operator()() const {return x;}
 };

 template<typename A> class const_qcache< A, 
  typename boost::enable_if<mpl::and_<is_matrix_expr<A>, mpl::not_<is_matrix_or_view<A> > > >::type> : boost::noncopyable { 
   typedef matrix<typename A::value_type> X;
   X x;
   public:
   const_qcache(A const & x_):x(x_){}
   typedef X const & exposed_type; 
   exposed_type operator()() const {return x;}
  };

 template<typename A> class const_qcache< A, 
  typename boost::enable_if<mpl::and_<is_vector_expr<A>, mpl::not_<is_vector_or_view<A> > > >::type> : boost::noncopyable { 
   typedef vector<typename A::value_type> X;
   X x;
   public:
   const_qcache(A const & x_):x(x_){}
   typedef typename X::view_type const & exposed_type; 
   exposed_type operator()() const {return x;}
  };

 template<typename A> class const_qcache< A, typename boost::enable_if<is_view_class<A> >::type> : boost::noncopyable { 
  typedef typename A::non_view_type data_type;
  boost::scoped_ptr< data_type  > data;
  boost::scoped_ptr< A > V_ptr;
  bool need_copy, init;
  A orig; A * V;
  void prep(){data.reset(new data_type(*V)); V_ptr.reset(new A(*data)); V=V_ptr.get(); init=true; } 
  public :
  explicit const_qcache(A const & x): need_copy (!( x.indexmap().is_contiguous())), init(false), orig(x), V(&orig){
   if (need_copy) std::cerr<< " I need a copy "<< std::endl;
  }
  const_qcache(const_qcache const &c):need_copy(c.need_copy),init(false),orig(c.orig),V(&orig) {}
  typedef A exposed_type; 
  exposed_type operator()() { if (need_copy && (!init)) prep(); return *(V);}
 };

 /*
  * reflexive_qcache.
  * Given A, a value or a view, it presents via the () operator
  *  - a & to the matrix if A is a value class (matrix, vector, array, ..)
  *  - a view if A is a view. If the view given at construction was not contiguous, it is a view to a 
  *    temporary regrouped value class. In that case, the data are back copied to the original at construction.
  *
  *  reflexive_qcache is NOT copyable. It should be used in local use only.
  *  [Rational : in simple case, like A =matrix, it allows to keep only a &, which is quicker
  *
  */
 template<typename A, typename Enable=void> struct reflexive_qcache;

 template<typename A> class reflexive_qcache< A, typename boost::enable_if<is_value_class<A> >::type> : boost::noncopyable { 
  A & x;
  public:
  reflexive_qcache(A & x_):x(x_){}
  typedef A & exposed_type; 
  exposed_type operator()() const {return x;}
 };

 template<typename A> class reflexive_qcache< A, typename boost::enable_if<is_view_class<A> >::type> : boost::noncopyable { 
  typedef typename A::non_view_type data_type;
  boost::scoped_ptr< data_type  > data;
  boost::scoped_ptr< A > V_ptr;
  bool need_copy, init;
  A orig; A * V;
  void prep(){data.reset(new data_type(*V)); V_ptr.reset(new A(*data)); V=V_ptr.get(); init=true; } 
  public :
  explicit reflexive_qcache(A const & x): need_copy (!( x.indexmap().is_contiguous())), init(false), orig(x), V(&orig){
   if (need_copy) std::cerr<< " I need a copy "<< std::endl;
  }
  reflexive_qcache(reflexive_qcache const &c):need_copy(c.need_copy),init(false),orig(c.orig),V(&orig) {}
  ~reflexive_qcache() { if (need_copy && init) orig = *(data);} // copy data back
  typedef A exposed_type; 
  exposed_type operator()() { if (need_copy && (!init)) prep(); return *(V);}
 };

}}//namespace triqs::arrays
#endif

