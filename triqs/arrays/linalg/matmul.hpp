
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

#ifndef TRIQS_ARRAYS_EXPRESSION_MATMUL_H
#define TRIQS_ARRAYS_EXPRESSION_MATMUL_H
#include "./matrix_cache.hpp"
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include <triqs/utility/proto/tools.hpp>

namespace triqs { namespace arrays { 
 namespace linalg {

  template<typename A, typename B> 
   class matmul_impl : 
    Tag::matrix_algebra_expression_terminal,Tag::expression_terminal, 
    Tag::has_special_assign, Tag::has_special_infix<'A'>, Tag::has_special_infix<'S'>, Tag::has_immutable_array_interface {

     //static_assert( (is_matrix_or_view<A>::value), "matmul : the first argument must be a matrix"); 
     //static_assert( (is_matrix_or_view<B>::value), "matmul : the second argument must be a matrix"); 
     typedef typename A::value_type V1;
     typedef typename B::value_type V2;
     static_assert((boost::is_same<V1,V2>::value),"Different values : not implemented");

     public:
     typedef BOOST_TYPEOF_TPL( V1() + V2()) value_type;
     typedef typename A::domain_type  domain_type;
     typedef typename domain_type::index_value_type index_value_type;
     A const & a; B const & b;
     // should be this for putting in expr... redo the cache....
     //typename A::view_type a; typename  B::view_type  b;
     //typename utility::proto::const_view_type_if_exists_else_type<A>::type a;
     //typename utility::proto::const_view_type_if_exists_else_type<B>::type b;
     mutable matrix_cache< A,matrix_cache_policy::const_copy, value_type> Ca;
     mutable matrix_cache< B,matrix_cache_policy::const_copy, value_type> Cb;

     private: 
     //typedef typename A::non_view_type matrix_type;
     //typedef typename matrix<value_type, typename A::opt_type> matrix_type;
     typedef matrix<value_type> matrix_type;
     mutable boost::shared_ptr<matrix_type> result;
     mutable bool init;
     void compute() const { 
      result.reset(new matrix_type(a.dim0(),  b.dim1()));
      std::cerr<<" -------------------\n matmul : compute"<<std::endl;
      std::cerr<<" a = "<< Ca() << "b = "<< Cb () <<std::endl;
      boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), 0.0, *result);
      std::cerr<<" a * b = "<< *result <<std::endl;
      init =true;
     }

     public:
     matmul_impl( A const & a_, B const & b_):a(a_),b(b_),Ca(a),Cb(b),init(false){
      if (a.dim1() != b.dim0()) TRIQS_RUNTIME_ERROR<< "Matrix product : dimension mismatch in A*B "<< a<<" "<< b; 
     }

     domain_type domain() const { return indexmaps::cuboid_domain<2>(mini_vector<size_t,2>(a.dim0(), b.dim1()));}
     size_t dim0() const { return a.dim0();} size_t dim1() const { return b.dim1();} 

     template<typename KeyType> value_type operator[] (KeyType const & key) const { if (!init) compute(); return (*result) [key]; }

     template<typename LHS> 
      void assign_invoke (LHS & lhs) const {// Optimized implementation of =
       static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix"); 
       std::cerr<<" -------------------\n matmul : assign_invoke"<<std::endl;
       std::cerr<<" a = "<< Ca() << "b = "<< Cb () <<std::endl;
       lhs.resize(dim0(),dim1());
       boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), 0.0, lhs);
       std::cerr<<" a * b = "<< lhs <<std::endl;
      }

     template<typename LHS> 
      void assign_add_invoke (LHS & lhs) const { //Optimized implementation of +=
       static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix"); 
       if (lhs.dim0() != dim0()) 
	TRIQS_RUNTIME_ERROR<< "Matmul : += operator : first dimension mismatch in A*B "<< lhs.dim0()<<" vs "<< dim0(); 
       if (lhs.dim1() != dim1()) 
	TRIQS_RUNTIME_ERROR<< "Matmul : += operator : first dimension mismatch in A*B "<< lhs.dim1()<<" vs "<< dim1(); 
       lhs.resize(dim0(),dim1());
       boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), 1.0, lhs);
      }

     template<typename LHS> 
      void assign_sub_invoke (LHS & lhs) const { //Optimized implementation of -=
       static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix"); 
       if (lhs.dim0() != dim0()) 
	TRIQS_RUNTIME_ERROR<< "Matmul : -= operator : first dimension mismatch in A*B "<< lhs.dim0()<<" vs "<< dim0(); 
       if (lhs.dim1() != dim1()) 
	TRIQS_RUNTIME_ERROR<< "Matmul : -= operator : first dimension mismatch in A*B "<< lhs.dim1()<<" vs "<< dim1(); 
       boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), -1.0, lhs);
      }

    };// class matmul_impl

  template<typename A, typename B> 
   std::ostream & operator<<(std::ostream & out, matmul_impl<A,B> const & x){ return out<<"matmul("<<x.a<<","<<x.b<<")";}

  template<typename A, typename B> matmul_impl<A,B> matmul (A const & a, B const & b) { return matmul_impl<A,B>(a,b); }

 } // linalg

 namespace result_of { template<typename A, typename B> struct matmul {  typedef linalg::matmul_impl<A,B> type;}; }

}}//namespace triqs::arrays 
#endif
