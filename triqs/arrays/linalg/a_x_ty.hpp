
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

#ifndef TRIQS_ARRAYS_EXPRESSION_A_X_TY_H
#define TRIQS_ARRAYS_EXPRESSION_A_X_TY_H
#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include "../matrix.hpp"
#include "../vector.hpp"
#include <boost/numeric/bindings/blas/level2/ger.hpp>

namespace triqs { namespace arrays { 

 namespace linalg { namespace details { template<typename ScalarType, typename VectorType1, typename VectorType2> class a_x_ty_impl; }}  // impl below
 namespace result_of { template<typename ScalarType,typename VectorType1, typename VectorType2> struct a_x_ty { 
  typedef linalg::details::a_x_ty_impl<ScalarType,VectorType1,VectorType2> type;}; 
 }

 namespace linalg { 

  template<typename ScalarType, typename VectorType1, typename VectorType2>
   details::a_x_ty_impl<ScalarType,VectorType1,VectorType2> a_x_ty (ScalarType a, VectorType1 const & x, VectorType2 const & y) { return details::a_x_ty_impl<ScalarType,VectorType1,VectorType2>(a,x,y); }

  //------------- IMPLEMENTATION -----------------------------------

  namespace details {

   template<typename ScalarType, typename VectorType1, typename VectorType2> 
    class a_x_ty_impl : TRIQS_MODEL_CONCEPT(ImmutableMatrix) { // TO BE DONE !!!! 
     // first check that VectorType1 and VectorType2 are matrices. 
     static_assert( (is_vector_or_view<VectorType1>::value), "a_x_ty : the first argument must be a vector"); 
     static_assert( (is_vector_or_view<VectorType2>::value), "a_x_ty : the second argument must be a vector"); 

     typedef typename VectorType1::value_type V1;
     typedef typename VectorType2::value_type V2;
     static_assert((boost::is_same<V1,V2>::value),"Different values : not implemented");

     public:
     typedef BOOST_TYPEOF_TPL( V1() + V2() + ScalarType()) value_type;
     typedef indexmaps::cuboid_domain<2> domain_type;
     ScalarType a;VectorType1 const & x; VectorType2 const & y; 

     private: 
     typedef typename VectorType2::non_view_type vector_type;

     public:
     a_x_ty_impl( ScalarType a_, VectorType1 const & x_, VectorType2 const & y_):a(a_),x(x_),y(y_){}
     domain_type domain() const { return domain_type(mini_vector<size_t,2>(x.dim0(), y.dim0()));}
     value_type operator[] (typename domain_type::index_value_type const & key) const { 
      return a * x(key[0]) * y(key[1]); }

     // Optimized implementation of =
     template<typename LHS> 
      friend void triqs_arrays_assign_delegation (LHS & lhs, a_x_ty_impl const & rhs)  {
       static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix or a matrix_view"); // check that the target is indeed a matrix.
       lhs()=0;
       boost::numeric::bindings::blas::ger(rhs.a, rhs.x, rhs.y,lhs);
      }

     //Optimized implementation of +=
     template<typename LHS> 
      friend void triqs_arrays_compound_assign_delegation (LHS & lhs, a_x_ty_impl const & rhs, mpl::char_<'A'>) {
       //std::cerr<<" Using optimized += for  a_x_ty_impl"<< std::endl;
       static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix or a matrix_view"); // check that the target is indeed a matrix.
       boost::numeric::bindings::blas::ger(rhs.a, rhs.x, rhs.y,lhs);
      }

     //Optimized implementation of -=
     template<typename LHS> 
      friend void triqs_arrays_compound_assign_delegation (LHS & lhs, a_x_ty_impl const & rhs, mpl::char_<'S'>) { 
       static_assert((is_matrix_or_view<LHS>::value), "LHS is not a matrix or a matrix_view"); // check that the target is indeed a matrix.
       boost::numeric::bindings::blas::ger(-rhs.a, rhs.x, rhs.y,lhs);
      }

     friend std::ostream & operator<<(std::ostream & out, a_x_ty_impl const & x){ return out<<"a_x_ty("<<x.a<<","<<x.x<<","<<x.y<<")";}

    };
  }
 } // linalg::details
}} // namespace triqs_arrays 
#endif
