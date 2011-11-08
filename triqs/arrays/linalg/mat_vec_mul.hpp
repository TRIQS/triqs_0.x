
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

#ifndef TRIQS_ARRAYS_EXPRESSION_MAT_VEC_MUL_H
#define TRIQS_ARRAYS_EXPRESSION_MAT_VEC_MUL_H

#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include "../matrix.hpp"
#include "../vector.hpp"
#include <boost/numeric/bindings/blas/level2/gemv.hpp>

namespace triqs { namespace arrays { 
 namespace linalg { namespace details { template<typename MT, typename VT> struct mat_vec_mul_impl; }}  // impl below
 namespace result_of { template<typename MT, typename VT> struct mat_vec_mul {  typedef linalg::details::mat_vec_mul_impl<MT,VT> type;}; }
}}

template<typename MT, typename VT>
typename boost::enable_if< boost::mpl::and_<triqs::arrays::is_matrix_or_view<MT>, triqs::arrays::is_vector_or_view<VT> >, 
	 triqs::arrays::linalg::details::mat_vec_mul_impl<MT,VT> >::type
operator* (MT const & a, VT const & b) { 
 return triqs::arrays::linalg::details::mat_vec_mul_impl<MT,VT>(a,b); 
}

namespace triqs { namespace arrays { namespace linalg {

 template<typename MT, typename VT> details::mat_vec_mul_impl<MT,VT> mat_vec_mul (MT const & a, VT const & b) { return details::mat_vec_mul_impl<MT,VT>(a,b); }

 namespace details { //------------- IMPLEMENTATION -----------------------------------

  template<typename MT, typename VT> 
   class mat_vec_mul_impl : Tag::expression_terminal, Tag::has_special_assign, Tag::has_special_infix<'A'>, Tag::has_special_infix<'S'>, Tag::has_immutable_array_interface {

    static_assert( (is_matrix_or_view<MT>::value), "mat_vec_mul : the first argument must be a matrix"); 
    static_assert( (is_vector_or_view<VT>::value), "mat_vec_mul : the second argument must be a vector"); 
    typedef typename MT::value_type V1;
    typedef typename VT::value_type V2;
    static_assert((boost::is_same<V1,V2>::value),"Different values : not implemented");

    public:
    typedef BOOST_TYPEOF_TPL( V1() + V2()) value_type;
    typedef typename VT::domain_type  domain_type;
    typedef typename domain_type::index_value_type index_value_type;
    MT const & M; VT const & V;

    private: 
    typedef typename VT::non_view_type vector_type;
    mutable boost::shared_ptr<vector_type> result;
    mutable bool init;
    void compute() const { 
     result.reset(new vector_type(M.dim0()));
     boost::numeric::bindings::blas::gemv(1, M, V, 0, *result);
     init =true;
    }

    public:
    mat_vec_mul_impl( MT const & M_, VT const & V_):M(M_),V(V_),init(false){}
    domain_type domain() const { return indexmaps::cuboid_domain<1>(mini_vector<size_t,1>(M.dim0()));}
    value_type operator[] (index_value_type const & key) const { if (!init) compute(); return (*result) [key]; }

    template<typename LHS> 
     void assign_invoke (LHS & lhs) const { // Optimized implementation of =
      static_assert((is_vector_or_view<LHS>::value), "LHS is not a vector or a vector_view"); 
      boost::numeric::bindings::blas::gemv(1, M, V, 0, lhs);
     }
    template<typename LHS> 
     void assign_add_invoke (LHS & lhs) const { // Optimized implementation of +=
      static_assert((is_vector_or_view<LHS>::value), "LHS is not a vector or a vector_view"); 
      boost::numeric::bindings::blas::gemv(1, M, V, 1, lhs);
     }
    template<typename LHS> 
     void assign_sub_invoke (LHS & lhs) const { //Optimized implementation of -=
      static_assert((is_vector_or_view<LHS>::value), "LHS is not a vector or a vector_view");
      boost::numeric::bindings::blas::gemv(1, M, V, -1, lhs);
     }
   };
  template<typename MT, typename VT> 
   std::ostream & operator<<(std::ostream & out, mat_vec_mul_impl<MT,VT> const & x){ return out<<"mat_vec_mul("<<x.M<<","<<x.V<<")";}
 }

}}}//namespace linalg::triqs::arrays 
#endif

