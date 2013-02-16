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
#ifndef TRIQS_ARRAYS_BLAS_LAPACK_GEMV_H
#define TRIQS_ARRAYS_BLAS_LAPACK_GEMV_H
#include "./tools.hpp"
#include "./qcache.hpp"

namespace triqs { namespace arrays { namespace blas { 

 using namespace blas_lapack_tools;
 namespace f77 { // overload
  typedef std::complex<double> dcomplex;

  extern "C" { 
   void TRIQS_FORTRAN_MANGLING(dgemv)(const char* trans, const int & m, const int & n, const double & alpha, const double A[], int & lda,
     const double x[], const int & incx, const double & beta, double y[], const int & incy);

   void TRIQS_FORTRAN_MANGLING(zgemv)(const char* trans, const int & m, const int & n, const std::complex<double> & alpha, const std::complex<double> A[], int & lda,
     const std::complex<double> x[], const int & incx, const std::complex<double> & beta, std::complex<double> y[], const int & incy);
  }

  void gemv (char * trans, const int & M, const int & N, double & alpha, const double* A, int & LDA, const double* x, const int & incx, double & beta, double* Y, const int & incy) { 
   TRIQS_FORTRAN_MANGLING(dgemv)(trans, M, N, alpha, A, LDA,x, incx , beta, Y, incy);
  }

  void gemv (char * trans, const int & M, const int & N, dcomplex & alpha, const dcomplex* A, int & LDA, const dcomplex* x, const int & incx, dcomplex & beta, dcomplex* Y, const int & incy) { 
   TRIQS_FORTRAN_MANGLING(zgemv)(trans, M, N, alpha, A, LDA,x, incx , beta, Y, incy);
  }
 }

 /**
  * Calls gemv
  * Takes care of making temporary copies if necessary
  */
 template<typename MT, typename VT, typename VTOut> 
  typename std::enable_if< is_blas_lapack_type<typename MT::value_type>::value && have_same_value_type< MT, VT, VTOut>::value >::type 
  gemv (typename MT::value_type alpha, MT const & A, VT const & X, typename MT::value_type beta, VTOut & Y) { 
   resize_or_check_if_view(Y,make_shape(A.dim0()));// first resize if necessary and possible 
   const_qcache<MT> Ca(A);
   const_qcache<VT> Cx(X); // mettre la condition a la main
   if (!(Ca().dim1() == Cx().size())) TRIQS_RUNTIME_ERROR << "Dimension mismatch in gemv : A : "<< Ca().shape() <<" while X : "<<Cx().shape();
   char trans_a= get_trans(Ca(), false); 
   int m1 = get_n_rows(Ca()), m2 = get_n_cols(Ca());
   int lda = get_ld(Ca());
   f77::gemv(&trans_a,m1,m2,alpha, Ca().data_start(), lda, Cx().data_start(), Cx().stride(), beta, Y.data_start(), Y.stride());
  }

}}}// namespace


#endif

