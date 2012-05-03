
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

#ifndef EXPRESSION_DETERMINANT_H
#define EXPRESSION_DETERMINANT_H

#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include "../matrix.hpp"
#include "../vector.hpp"
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>

namespace triqs { namespace arrays { 
 namespace linalg { namespace details { template<typename A> struct determinant_impl; }} // impl below

 namespace result_of { 
  template<typename A> struct determinant {  
   static_assert( (is_matrix_or_view<A>::value), "determinant : the argument must be a matrix"); 
   typedef linalg::details::determinant_impl<A> type;}; 
 }

 namespace linalg {

  template<typename A> 
   typename result_of::determinant<A>::type determinant (A const & a) { 
    return details::determinant_impl<A>(a); 
   }

  //------------- IMPLEMENTATION -----------------------------------

  namespace details {

   //using namespace indexmaps;

   template<typename A> 
    class determinant_impl : Tag::expression_terminal, Tag::has_immutable_array_interface {
     typedef typename A::non_view_type matrix_type;
     public:
     typedef typename A::value_type value_type;
     typedef indexmaps::cuboid_domain<0> domain_type;
     A const & a; 

     static value_type compute_det_after_getrf(matrix_type const & M,  triqs::arrays::vector <int> const & ipiv) {
      const size_t n = M.dim0();
      value_type det =1;
      for (size_t i =0; i<n; i++) det *= M(i,i);
      bool flip=false;// compute the sign of the permutation
      for (size_t i=0; i<n; i++) {if (ipiv(i)!=int(i)+1) flip = !(flip);}
      return (flip ? - det : det) ;
     }

     private: 
     typedef typename domain_type::index_value_type index_value_type;
     mutable boost::shared_ptr<matrix_type> m_copy;
     mutable value_type result;
     mutable bool init;

     static value_type invoke( matrix_type & M)  {
      const size_t n = M.dim0();
      if (n==0) return 1;
      triqs::arrays::vector <int> ipiv(n);
      int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
      TRIQS_ARRAYS_CHECK_OR_THROW(( info==0), "Matrix inversion impossible for"<<M<<"lapack has returned error :"<<info);
      return compute_det_after_getrf(M,ipiv);
     }

     void compute() const { 
      m_copy.reset(new matrix_type(a)); // makes a copy
      result = invoke(*m_copy); //do it
      init =true;
     }

     public:
     determinant_impl( A const & a_):a(a_),init(false){} //Check::is_square_matrix(a,"determinant");}
     domain_type domain() const { return domain_type(); } 
     value_type operator[] (index_value_type const & key) const { if (!init) compute(); return result; }
     operator value_type() const { if (!init) compute(); return result; };
    };

   template<typename A> 
    std::ostream & operator<<(std::ostream & out, determinant_impl<A> const & x){ return out<<"determinant("<<x.a<<")";}
  } // details
 } // linalg
}} // namespace triqs_arrays 
#endif

