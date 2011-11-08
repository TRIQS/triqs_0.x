
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

#ifndef EXPRESSION_INVERSE_H
#define EXPRESSION_INVERSE_H

#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/enable_if.hpp>
#include "../impl/common.hpp"
#include "../matrix.hpp"
#include "../vector.hpp"
#include "./determinant.hpp"
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>

namespace triqs { namespace arrays { 
 namespace linalg { namespace details { template<typename A,bool compute_determinant> struct inverse_impl; }}  // impl below

 namespace result_of { 
  template<typename A> struct inverse {  
   static_assert( (is_matrix_or_view<A>::value), "inverse : the argument must be a matrix"); 
   typedef linalg::details::inverse_impl<A,false> type;
  }; 
  template<typename A> struct inverse_and_compute_det {
   static_assert( (is_matrix_or_view<A>::value), "inverse_and_compute_det : the argument must be a matrix"); 
   typedef linalg::details::inverse_impl<A,true> type;}; 
 }

 namespace linalg { 
  template<typename A> details::inverse_impl<A,false> inverse (A const & a) { return details::inverse_impl<A,false>(a); }

  template<typename A> details::inverse_impl<A,true> inverse_and_compute_det (A const & a) { return details::inverse_impl<A,true>(a); }

  //------------- IMPLEMENTATION -----------------------------------

  namespace details {

   // this class is to add the possibility of computing determinant.
   // using CRTP as a "static" virtual function
   template<typename Derived, typename A, bool CD> struct inverse_impl_det_access {
    void compute_det(typename A::non_view_type const & M,  triqs::arrays::vector <int> const & ipiv) const {}
   };

   template<typename Derived, typename A> struct inverse_impl_det_access<Derived,A,true> {
    mutable typename A::value_type det;
    void compute_det(typename A::non_view_type const & M,  triqs::arrays::vector <int> const & ipiv) const {
     det = determinant_impl<A>::compute_det_after_getrf(M,ipiv);
    }
    typename A::value_type determinant() const { 
     Derived const & self = static_cast<Derived const & >(*this); 
     if (!self.init) self.compute(); 
     return det; 
    }
   };

   //-----------------------------------------------------------

   template<typename A,bool compute_determinant>
    class inverse_impl : Tag::expression_terminal, Tag::has_special_assign, Tag::has_immutable_array_interface,
     public inverse_impl_det_access<inverse_impl<A,compute_determinant>,A,compute_determinant>
   {
    public:
     typedef typename A::value_type value_type;
     typedef typename A::domain_type domain_type;
     typedef typename domain_type::index_value_type index_value_type;
     A const & a; 

     inverse_impl( A const & a_):a(a_),init(false){ }// Check::is_square_matrix(a,"inverse");}
     domain_type domain() const { return a.domain(); } 
     value_type operator[] (index_value_type const & key) const { if (!init) compute(); return (*result) [key]; }

     template<typename Opt> 
      void assign_invoke (matrix<value_type,Opt> & lhs) const { // Optimized implementation of =
       if ((matrix<value_type,Opt>::order !=a.order) || (lhs.data_start() != a.data_start())) {lhs = a;} 
       // if M = inverse(M) with the SAME object, then we do not need to copy the data
       invoke(lhs);
      }
     template<typename Opt> 
      void assign_invoke (matrix_view<value_type,Opt> & lhs) const { // Optimized implementation of =
       if ((matrix<value_type,Opt>::order  !=a.order)|| (lhs.data_start() != a.data_start())) {lhs = a;} 
       // if M = inverse(M) with the SAME object, then we do not need to copy the data
       invoke(lhs);
      }

    private:
     friend class inverse_impl_det_access<inverse_impl<A,compute_determinant>,A,compute_determinant>;
     typedef typename A::non_view_type matrix_type;
     mutable boost::shared_ptr<matrix_type> result;
     mutable bool init;
     template<typename Opt> 
     void invoke(matrix<value_type,Opt> & M) const { matrix_view<value_type,Opt> V(M);invoke(V);} 
     template<typename ST> 
      void invoke(matrix_view<value_type,Option::options<Tag::C,ST> > & M) const {matrix_view<value_type, Option::options<Tag::Fortran,ST> > V(M.transpose());invoke(V);} 
     template<typename ST> 
      void invoke(matrix_view<value_type,Option::options<Tag::Fortran,ST> > & M) const { 
       triqs::arrays::vector <int> ipiv(M.dim0());
       int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
       if (info!=0) TRIQS_ARRAYS_THROW("Matrix "<<M<<"is not invertible");
       this->compute_det(M,ipiv);
       info = boost::numeric::bindings::lapack::getri(M, ipiv);
      }
     void compute() const { 
      result.reset(new matrix_type(a)); // makes a copy
      invoke(*result); //do it
      init =true;
     }
   };

   template<typename A, bool CD> 
    std::ostream & operator<<(std::ostream & out, inverse_impl<A,CD> const & x){ return out<<"inverse("<<x.a<<")";}
  }} // linalg::details
}} // namespace triqs_arrays 

#endif

