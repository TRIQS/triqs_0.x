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
#ifndef TRIQS_ARRAYS_LINALG_DET_INV_H
#define TRIQS_ARRAYS_LINALG_DET_INV_H

#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/enable_if.hpp>
#include "../impl/common.hpp"
#include "../qcache.hpp"
#include "triqs/utility/proto/tools.hpp"
#include "../matrix.hpp"
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>

namespace triqs { namespace arrays { namespace linalg { 

 /// Error which occurs during the matrix inversion 
 class matrix_inverse_exception : public triqs::runtime_error {};

 /**
  * Lazy result of inverse(M) where M can be : 
  *  * a matrix, a matrix_view
  *  * any matrix expression
  * The object is lazy, it does not do any computation.
  * It can be copied at no cost
  * It keeps view of the object A if it a matrix, a copy if it is a formal expression.
  */ 
 template<typename A, class Enable = void> class inverse_lazy;

 /// Lazy inversion 
 template<class A> inverse_lazy<A> inverse (A const & a) { return inverse_lazy<A>(a); }
 
 ///
 template<typename A> struct determinant_lazy;
 
 /// Lazy computation of det
 template<typename A> determinant_lazy<A> determinant (A const & a) { return determinant_lazy<A>(a); }

 // ----------------- implementation -----------------------------------------

 //worker takes a contiguous view and compute the det and inverse in two steps.
 //it is separated in case of multiple use (no reallocation of ipvi, etc...)
 //NB a view does not resize, only its elements can be changed
 template<typename ViewType> class det_and_inverse_worker { 
  static_assert ( (is_matrix_view<ViewType>::value),"class must have be a view");
  typedef typename ViewType::value_type VT;
  typedef matrix_view<VT,Option::Fortran > V_type;
  ViewType V;
  const size_t dim; 
  triqs::arrays::vector <int> ipiv;
  short step; 

  public:
  det_and_inverse_worker (ViewType const & a): V(a), dim(a.dim0()), ipiv(dim), step(0) { 
   if (a.dim0()!=a.dim1()) TRIQS_RUNTIME_ERROR<<"Inverse/Det error : non-square matrix. Dimensions are : ("<<a.dim0()<<","<<a.dim1()<<")"<<"\n  ";
   if (!(has_contiguous_data(a))) TRIQS_RUNTIME_ERROR<<"det_and_inverse_worker only takes a contiguous view";
  }
  VT det() { V_type W = fortran_view(V); _step1(W); _compute_det(W); return _det;}
  ViewType const & inverse() { if (step<2) { V_type W = fortran_view(V); _step1(W); _step2(W);} return V;}

  private:
  int info; VT _det;

  matrix_view<VT,Option::Fortran> fortran_view (matrix<VT,Option::C> const & x) {return x.transpose();}
  matrix_view<VT,Option::Fortran> fortran_view (matrix<VT,Option::Fortran> const & x) {return x;}
  matrix_view<VT,Option::Fortran> fortran_view (matrix_view<VT,Option::C> const & x) {return x.transpose();}
  matrix_view<VT,Option::Fortran> fortran_view (matrix_view<VT,Option::Fortran> const & x) {return x;}

  void _step1(V_type & W) { 
   if (step >0) return;
   step=1; 
   info = boost::numeric::bindings::lapack::getrf(W, ipiv);
   if (info<0) throw matrix_inverse_exception() << "Inverse/Det error : failure of getrf lapack routine ";
  }

  void _compute_det(V_type const & W) { 
   if (step>1) return;
   _det =1;
   for (size_t i =0; i<dim; i++) _det *= W(i,i);
   bool flip=false;// compute the sign of the permutation
   for (size_t i=0; i<dim; i++) {if (ipiv(i)!=int(i)+1) flip = !(flip);}
   _det= (flip ? - _det : _det) ;
  }

  void _step2(V_type & W) { 
   assert(step==1); //if (step==1) return;
   step=2;
   _compute_det(W); 
   info = boost::numeric::bindings::lapack::getri(W, ipiv);
   std::cerr<<" after tri "<< W<<std::endl;
   if (info!=0) throw matrix_inverse_exception() << "Inverse/Det error : matrix is not invertible";
  }
 };

 //-----------------------------------------------------------

 // an implementation class to gather the common part to matrix and expression....
 template<typename A> struct inverse_lazy_impl : Tag::matrix_algebra_expression_terminal, Tag::has_immutable_array_interface {
  typedef typename A::value_type value_type;
  typedef typename A::domain_type domain_type;
  typedef typename domain_type::index_value_type index_value_type;
  typedef typename utility::proto::const_view_type_if_exists_else_type<A>::type A_type;
  const A_type a;
  inverse_lazy_impl(A const & a_):a (a_)  {}
  domain_type domain() const { return a.domain(); } 
  value_type operator[] (index_value_type const & key) const { activate();  _id->M [key]; }
  friend std::ostream & operator<<(std::ostream & out,inverse_lazy_impl const&x){return out<<"inverse("<<x.a<<")";}
  protected: 
  struct internal_data { // implementing the pattern LazyPreCompute
   typedef typename A_type::non_view_type M_type;
   M_type M;
   internal_data(inverse_lazy_impl const & P):M(P.a){det_and_inverse_worker<A_type> worker(M); worker.inverse();}
  };
  friend struct internal_data;
  mutable boost::shared_ptr<internal_data> _id;
  void activate() const { if (!_id) _id= boost::make_shared<internal_data>(*this);}
 };

 // The general case
 template<typename A>
  struct inverse_lazy<A,typename boost::disable_if< is_matrix_or_view<A> >::type > : inverse_lazy_impl<A> { 
   inverse_lazy(A const & a):inverse_lazy_impl<A>(a) { }
  };

 // for matrix and matrix_views, we have more optimisation possible ....
 template<typename A>
  struct inverse_lazy<A,typename boost::enable_if< is_matrix_or_view<A> >::type >:inverse_lazy_impl<A>, Tag::has_special_assign{
   inverse_lazy(A const & a):inverse_lazy_impl<A>(a) { }

   template<typename MT> 
    void assign_invoke (MT & lhs) const { // Optimized implementation of =
     static_assert(is_matrix_or_view<MT>::value, "Internal error");
     if ((MT::order  !=this->a.order)|| (lhs.data_start() != this->a.data_start()) || !(has_contiguous_data(lhs))) { this->activate(); lhs = this->_id->M;} 
     else {// if M = inverse(M) with the SAME object, then we do not need to copy the data
      reflexive_qcache<MT> C(lhs);// a reflexive cache will use a temporary "regrouping" copy if and only if needed
      det_and_inverse_worker<typename MT::view_type> W(C());// the worker to make the inversion of the lhs... 
      W.inverse(); // worker is working ...
     }
    }
   friend std::ostream & operator<<(std::ostream & out,inverse_lazy const&x){return out<<"inverse("<<x.a<<")";}
  };

 //------------------- det   ----------------------------------------

 template<typename A> struct determinant_lazy : Tag::expression_terminal, Tag::scalar_expression_terminal {
   typedef typename A::value_type value_type;
   typedef typename utility::proto::const_view_type_if_exists_else_type<A>::type A_type;
   A_type a;
   determinant_lazy(A const & a_):a(a_){}
   operator value_type()  { activate(); return _id->det; }
   friend std::ostream & operator<<(std::ostream & out, determinant_lazy const & x){ return out<<"determinant("<<x.a<<")";}
   protected:
   struct internal_data {
    typedef typename A_type::non_view_type M_type;
    M_type M; typename A::value_type det;
    internal_data(determinant_lazy const & P):M(P.a){det_and_inverse_worker<A_type> worker(M); det = worker.det();}
   };
   friend struct internal_data;
   mutable boost::shared_ptr<internal_data> _id;
   void activate() const { if (!_id) _id= boost::make_shared<internal_data>(*this);}
  };


}}} // namespace triqs::arrays::linalg
#endif
