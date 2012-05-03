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
#include "../qcache.hpp"
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include <triqs/utility/proto/tools.hpp>

namespace triqs { namespace arrays { namespace linalg {

 template<typename A, typename B> 
  class matmul_lazy : 
   Tag::matrix_algebra_expression_terminal,Tag::expression_terminal, 
   Tag::has_special_assign, Tag::has_special_infix<'A'>, Tag::has_special_infix<'S'>, Tag::has_immutable_array_interface {

    typedef typename A::value_type V1;
    typedef typename B::value_type V2;
    static_assert((boost::is_same<V1,V2>::value),"Different values : not implemented");

    public:
    typedef BOOST_TYPEOF_TPL( V1() * V2()) value_type;
    typedef typename A::domain_type  domain_type;
    typedef typename utility::proto::const_view_type_if_exists_else_type<A>::type A_type;
    typedef typename utility::proto::const_view_type_if_exists_else_type<B>::type B_type;
    const A_type a; const B_type b;

    private: 
    typedef matrix<value_type> matrix_type;

    struct internal_data {
     matrix_type R;
     internal_data(matmul_lazy const & P): R( P.a.dim0(), P.b.dim1()) {
      const_qcache<A_type> Ca(P.a); const_qcache<B_type> Cb(P.b);
      std::cerr<<" -------------------\n matmul : compute"<<std::endl;
      std::cerr<<" a = "<< Ca() << "b = "<< Cb () <<std::endl;
      boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), 0.0, R);
      std::cerr<<" a * b = "<< R <<std::endl;
     }
    };
    friend struct internal_data;
    mutable boost::shared_ptr<internal_data> _ip;
    void activate() const { if (!_ip) _ip= boost::make_shared<internal_data>(*this);}

    public:
    matmul_lazy( A const & a_, B const & b_):a(a_),b(b_){ 
     if (a.dim1() != b.dim0()) TRIQS_RUNTIME_ERROR<< "Matrix product : dimension mismatch in A*B "<< a<<" "<< b; 
    }

    domain_type domain() const { return indexmaps::cuboid_domain<2>(mini_vector<size_t,2>(a.dim0(), b.dim1()));}
    size_t dim0() const { return a.dim0();} 
    size_t dim1() const { return b.dim1();} 

    template<typename KeyType> value_type operator[] (KeyType const & key) const { activate(); return _ip->R [key]; }

    template<typename LHS> 
     void assign_invoke (LHS & lhs) const {// Optimized implementation of =
      static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix"); 
      const_qcache<A_type> Ca(a); const_qcache<B_type> Cb(b);
      std::cerr<<" -------------------\n matmul : assign_invoke"<<std::endl;
      std::cerr<<" a = "<< Ca() << "b = "<< Cb () <<std::endl;
      resize_or_check_if_view(lhs,make_shape(dim0(),dim1()));
      //lhs.resize(dim0(),dim1());
      boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), 0.0, lhs);
      std::cerr<<" a * b = "<< lhs <<std::endl;
     }

    template<typename LHS> void assign_add_invoke (LHS & lhs) const { assign_comp_impl(lhs,1.0);}
    template<typename LHS> void assign_sub_invoke (LHS & lhs) const { assign_comp_impl(lhs,-1.0);}

     private:   
    template<typename LHS> void assign_comp_impl (LHS & lhs, double S) const { 
     static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix"); 
     if (lhs.dim0() != dim0()) 
      TRIQS_RUNTIME_ERROR<< "Matmul : -= operator : first dimension mismatch in A*B "<< lhs.dim0()<<" vs "<< dim0(); 
     if (lhs.dim1() != dim1()) 
      TRIQS_RUNTIME_ERROR<< "Matmul : -= operator : first dimension mismatch in A*B "<< lhs.dim1()<<" vs "<< dim1(); 
     const_qcache<A_type> Ca(a); const_qcache<B_type> Cb(b);
     boost::numeric::bindings::blas::gemm(1.0,Ca(), Cb(), S, lhs);
    }

    friend std::ostream & operator<<(std::ostream & out, matmul_lazy<A,B> const & x){return out<<x.a<<" * "<<x.b;}
    //friend std::ostream & operator<<(std::ostream & out, matmul_lazy<A,B> const & x){return out<<"matmul("<<x.a<<","<<x.b<<")";}

    };// class matmul_lazy

    template<typename A, typename B> matmul_lazy<A,B> matmul (A const & a, B const & b) { return matmul_lazy<A,B>(a,b); }

   }}}//namespace triqs::arrays::linalg
#endif
