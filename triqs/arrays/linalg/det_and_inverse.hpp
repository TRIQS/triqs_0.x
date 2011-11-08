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
#include "./matrix_cache.hpp"
//#include "../vector.hpp"
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>

namespace triqs { namespace arrays { namespace linalg { 

 class matrix_inverse_exception : public triqs::runtime_error {};

  template<class T> matrix_view<T,Option::Fortran> fortran_view (matrix<T,Option::C> const & x) {return x.transpose();}
  template<class T> matrix_view<T,Option::Fortran> fortran_view (matrix<T,Option::Fortran> const & x) {return x;}
  template<class T> matrix_view<T,Option::Fortran> fortran_view (matrix_view<T,Option::C> const & x) {return x.transpose();}
  template<class T> matrix_view<T,Option::Fortran> fortran_view (matrix_view<T,Option::Fortran> const & x) {return x;}

  template<typename A, matrix_cache_policy::policy_type cache_policy>
   class det_and_inverse_worker { 
    static_assert ( (has_immutable_array_interface<A>::value),"Class must have ImmutableArrayInterface concept");
    typedef typename A::value_type value_type;
    typedef matrix_cache< A, cache_policy, value_type> cache_type;
    typedef typename cache_type::exposed_type exposed_type ;
    typedef matrix_view<value_type,Option::Fortran > V_type;
    cache_type C; 
    const size_t dim; 
    triqs::arrays::vector <int> ipiv;
    short step; 

    public:
    det_and_inverse_worker (typename cache_type::constructor_arg_type a): C(a), dim(a.dim0()), ipiv(dim), step(0) { 
     if (a.dim0()!=a.dim1()) TRIQS_RUNTIME_ERROR<<"Inverse/Det error : non-square matrix. Dimensions are : ("<<a.dim0()<<","<<a.dim1()<<")"<<"\n  ";
    }
    value_type det() { V_type W = fortran_view(C()); _step1(W); _compute_det(W); return _det;}
    exposed_type inverse() { if (step<2) { V_type W = fortran_view(C()); _step1(W); _step2(W);} return C();}

    private:
    int info; value_type _det;

    void _step1(V_type & W) { 
     if (step >0) return;
     step=1; 
     info = boost::numeric::bindings::lapack::getrf(W, ipiv);
     std::cerr<<" after tf "<< W<<std::endl;
     if (info<0) throw matrix_inverse_exception() << "Inverse/Det error : failure of getrf lapack routine ";
    }

    void _compute_det(V_type const & W) { 
     if (step>1) return;
     _det =1;
     for (size_t i =0; i<dim; i++) _det *= W(i,i);
     bool flip=false;// compute the sign of the permutation
     for (size_t i=0; i<dim; i++) {if (ipiv(i)!=int(i)+1) flip = !(flip);}
     _det= (flip ? - _det : _det) ;
     std::cerr<<" det is  "<< _det <<std::endl;
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

  template<typename A, class Enable = void> class inverse_expr;

  // for general expression.
  template<typename A>
   class inverse_expr<A,typename boost::disable_if< is_matrix_or_view<A> >::type > :
   Tag::matrix_algebra_expression_terminal, Tag::has_immutable_array_interface
   {
    det_and_inverse_worker<A,matrix_cache_policy::copy> worker;
    public:
    typedef typename A::value_type value_type;
    typedef typename A::domain_type domain_type;
    typedef typename domain_type::index_value_type index_value_type;
    A const & a; 

    inverse_expr( A const & a_):worker(a_),a (a_)  {}
    domain_type domain() const { return a.domain(); } 
    value_type operator[] (index_value_type const & key) const { return worker.inverse() [key]; }
   };

  // for matrix and matrix_views, we have more optimisation possible ....
  template<typename A>
   class inverse_expr<A,typename boost::enable_if< is_matrix_or_view<A> >::type > :
   Tag::matrix_algebra_expression_terminal, Tag::has_special_assign, Tag::has_immutable_array_interface {
    mutable det_and_inverse_worker<A,matrix_cache_policy::copy> worker;
    public:
    typedef typename A::value_type value_type;
    typedef typename A::domain_type domain_type;
    typedef typename domain_type::index_value_type index_value_type;
    A const & a; 

    inverse_expr( A const & a_):worker(a_),a (a_) {}
    domain_type domain() const { return a.domain(); } 
    value_type operator[] (index_value_type const & key) const { return worker.inverse() [key]; }

    template<typename Opt> 
     void assign_invoke (matrix<value_type,Opt> & lhs) const { // Optimized implementation of =
      if ((matrix<value_type,Opt>::order !=a.order) || (lhs.data_start() != a.data_start())) {lhs = worker.inverse();} 
      else {// if M = inverse(M) with the SAME object, then we do not need to copy the data 
       det_and_inverse_worker<matrix<value_type,Opt>, matrix_cache_policy::reflexive> W(lhs); 
       W.inverse();
      }
     }
    template<typename Opt> 
     void assign_invoke (matrix_view<value_type,Opt> & lhs) const { // Optimized implementation of =
      if ((matrix<value_type,Opt>::order  !=a.order)|| (lhs.data_start() != a.data_start())
	|| !(lhs.indexmap().is_contiguous())) {lhs = worker.inverse() ;} 
      else { // if M = inverse(M) with the SAME object, and a contiguous data, then we do not need to copy the data
       det_and_inverse_worker<matrix<value_type,Opt>, matrix_cache_policy::reflexive> W(lhs); 
       W.inverse();
      }
     }
   };

  template<class A> std::ostream & operator<<(std::ostream & out,inverse_expr<A> const&x){return out<<"inverse("<<x.a<<")";}
  template<class A> inverse_expr<A> inverse (A const & a) { return inverse_expr<A>(a); }
  namespace result_of { template<class A> struct inverse {  typedef inverse_expr<A>  type; }; }

  //-----------------------------------------------------------

  template<typename A> 
   class determinant_expr : Tag::expression_terminal, Tag::scalar_expression_terminal {
    det_and_inverse_worker<A,matrix_cache_policy::copy> worker;
    public:
    typedef typename A::value_type value_type;
    A const & a; 
    determinant_expr(A const & a_):worker(a_),a(a_){}
    operator value_type()  { return worker.det(); }
   };

  template<typename A> std::ostream & operator<<(std::ostream & out, determinant_expr<A> const & x){ return out<<"determinant("<<x.a<<")";}
  template<typename A> determinant_expr<A> determinant (A const & a) { return determinant_expr<A>(a); }
  namespace result_of { template<typename A> struct determinant {  typedef determinant_expr<A>  type; }; }

}}} // namespace triqs::arrays::linalg
#endif
