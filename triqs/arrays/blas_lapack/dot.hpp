/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_BLAS_LAPACK_DOT_H
#define TRIQS_ARRAYS_BLAS_LAPACK_DOT_H
#include <triqs/utility/fortran_mangling.hpp>
#include "./is_blas_lapack_type.hpp"
//#include "./qcache.hpp"

namespace triqs { namespace arrays { namespace blas { 

 namespace f77 { // overload

  extern "C" { 
   double TRIQS_FORTRAN_MANGLING(ddot)(const int & N, const double * x, const int& incx, const double * y, const int & incy);
   double TRIQS_FORTRAN_MANGLING(zdotc)(const int & N, const std::complex<double> * x, const int& incx, const std::complex<double> * y, const int& incy);
  }

  double dot (const int & M, const double* x, const int & incx, const double* Y, const int & incy)  { 
   return TRIQS_FORTRAN_MANGLING(ddot)(M, x, incx,  Y, incy);
  }
  std::complex<double> dot (const int & M, const std::complex<double>* x, const int & incx, const std::complex<double>* Y, const int & incy)  { 
   return TRIQS_FORTRAN_MANGLING(zdotc)(M, x, incx,  Y, incy);
  }
 }

 /**
  * Calls dot product of 2 vectors.
  * Takes care of making temporary copies if necessary
  */
 template< typename VectorXType,  typename VectorYType> 
  typename std::enable_if< is_blas_lapack_type<typename VectorXType::value_type>::value && have_same_value_type< VectorXType, VectorYType>::value, double >::type 
  dot (VectorXType const & X, VectorYType const & Y) { 
   static_assert( is_amv_value_or_view_class<VectorXType>::value, "blas1 bindings only take vector and vector_view");
   static_assert( is_amv_value_or_view_class<VectorYType>::value, "blas1 bindings only take vector and vector_view");
   if (( X.size() != Y.size()) ) TRIQS_RUNTIME_ERROR << "Dimension mismatch in dot : X : "<<X().shape()<<" and Y : "<<Y().shape();
   //const_qcache<VectorXType> Cx(X); // mettre la condition a la main
   //const_qcache<VectorYType> Cy(Y); // mettre la condition a la main
   return f77::dot(X.size(), X.data_start(), X.stride(), Y.data_start(), Y.stride());
   //return f77::dot(X.size(), Cx().data_start(), Cx().stride(), Cy().data_start(), Cy().stride());
  }

}}}// namespace

#endif

