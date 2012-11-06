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
namespace triqs { namespace arrays { 

 /**
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
 template<typename A, typename Enable=void> class const_qcache;

 template<typename A> class const_qcache< A, ENABLE_IF(is_amv_value_class<A>) > : boost::noncopyable { 
   A const & x;
   public:
   const_qcache(A const & x_):x(x_){}
   typedef A const & exposed_type; 
   exposed_type operator()() const {return x;}
  };

 template<typename A> class const_qcache< A, ENABLE_IF(mpl::and_<ImmutableMatrix<A>, mpl::not_<is_matrix_or_view<A> > >) > : 
  boost::noncopyable { 
   typedef matrix<typename A::value_type> X;
   X x;
   public:
   const_qcache(A const & x_):x(x_){}
   typedef X const & exposed_type; 
   exposed_type operator()() const {return x;}
  };

 template<typename A> class const_qcache< A, ENABLE_IF(mpl::and_<ImmutableVector<A>, mpl::not_<is_vector_or_view<A> > >) > :
  boost::noncopyable { 
   typedef vector<typename A::value_type> X;
   X x;
   public:
   const_qcache(A const & x_):x(x_){}
   typedef X const & exposed_type; 
   exposed_type operator()() const {return x;}
  };

 template<typename A> class const_qcache< A, ENABLE_IF(is_amv_view_class<A>)> : boost::noncopyable { 
  const bool need_copy;
  A keeper;
  struct internal_data {
   typename A::non_view_type copy_data;
   A view;
   internal_data(const_qcache const & P) : copy_data(P.keeper), view(copy_data) {
#ifdef TRIQS_ARRAYS_CACHE_COPY_VERBOSE
    std::cerr<< " Cache : copy made "<< std::endl;
    std::cerr<< " -- TRACE = --" << std::endl << triqs::utility::stack_trace() << std::endl;
#endif
   }
  };
  friend struct internal_data;   
  mutable boost::scoped_ptr<internal_data> _id;   
  internal_data & id() const { assert(need_copy); if (!_id) _id.reset(new internal_data(*this)); return *_id;}

  protected:
  void back_update() { if (_id) keeper = _id->copy_data; } // copy data back (for derivative, but but mechanism here).
  A & view() { return (need_copy ? id().view : keeper);}
  A const & view() const{ return (need_copy ? id().view : keeper);}

  public :
  explicit const_qcache(A const & x): need_copy (!(has_contiguous_data(x))), keeper(x) {}
  typedef const A exposed_type; 
  exposed_type operator()() const { return view();}
 };

 /**
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
 template<typename A, typename Enable=void> class reflexive_qcache;

 template<typename A> class reflexive_qcache <A, ENABLE_IF(is_amv_value_class<A>)> : boost::noncopyable { 
  A & x;
  public:
  reflexive_qcache(A & x_):x(x_){}
  typedef A & exposed_type; 
  exposed_type operator()() const {return x;}
 };

 template<typename A> class reflexive_qcache <A, ENABLE_IF(is_amv_view_class<A>)> : const_qcache<A> {
  typedef const_qcache<A> B;
  public :
  explicit reflexive_qcache(A const & x) : B(x) {} 
  ~reflexive_qcache() { this->back_update();}
  typedef A exposed_type; 
  exposed_type operator()() { return B::view(); } 
 };

}}//namespace triqs::arrays
#endif

