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
#ifndef TRIQS_ARRAYS_BLAS_LAPACK_GEMM_H
#define TRIQS_ARRAYS_BLAS_LAPACK_GEMM_H
#include <complex>
#include <triqs/utility/fortran_mangling.hpp>
#include "./is_blas_lapack_type.hpp"
#include "./qcache.hpp"

namespace triqs { namespace arrays { namespace blas { 
 
 using namespace blas_lapack_tools;
 namespace f77 { // overload

  extern "C" { 
   void TRIQS_FORTRAN_MANGLING(dgemm) (char *, char *,  int & , int & , int & , double &, const double[], int &, const double[], int &, double &, double[], int & ); 
   void TRIQS_FORTRAN_MANGLING(zgemm) (char *, char *,  int & , int & , int & , std::complex<double> &, const std::complex<double>[], int &,
     const std::complex<double>[], int &, std::complex<double> &, std::complex<double>[], int & );
  }

  // first version of gemm : just to fuse dgemm and zgemm ...
  void gemm (char * trans_a, char *trans_b, int & M, int & N, int & K, double & alpha, const double* A, int & LDA, const double* B, int & LDB, double & beta, double* C, int &  LDC) { 
   TRIQS_FORTRAN_MANGLING(dgemm)(trans_a,trans_b,M,N,K,alpha, A, LDA, B, LDB, beta, C, LDC);
  }

  typedef std::complex<double> dcomplex;
  void gemm (char * trans_a, char *trans_b, int & M, int & N, int & K, dcomplex & alpha, const dcomplex* A, int & LDA, const dcomplex* B, int & LDB, dcomplex & beta, dcomplex* C, int &  LDC) { 
   TRIQS_FORTRAN_MANGLING(zgemm)(trans_a,trans_b,M,N,K,alpha, A, LDA, B, LDB, beta, C, LDC);
  }
 }

 /**
  * Calls gemm on a matrix or view
  * Takes care of making temporary copies if necessary
  */
 template<typename MatrixType1, typename MatrixType2, typename MatrixTypeOut> 
  typename std::enable_if< is_blas_lapack_type<typename MatrixType1::value_type>::value && have_same_value_type< MatrixType1, MatrixType2, MatrixTypeOut>::value >::type 
  gemm (typename MatrixType1::value_type alpha, MatrixType1 const & A, MatrixType2 const & B, typename MatrixType1::value_type beta, MatrixTypeOut & C) { 

   // first resize if necessary and possible 
   resize_or_check_if_view(C,make_shape(A.dim0(),B.dim1()));

   // now we use qcache instead of the matrix to make a copy if necessary ... 
   // not optimal : if stride == 1, N ---> use LDA parameters
   // change the condition in the qcache construction....
   reflexive_qcache<MatrixTypeOut> Cc(C);

   if (C.memory_layout_is_c()) {
    // then tC = tB tA ! 
    const_qcache<MatrixType1> Cb(A); // note the inversion  A <-> B
    const_qcache<MatrixType2> Ca(B); // note the inversion  A <-> B
    if (!(Ca().dim0() == Cb().dim1())) TRIQS_RUNTIME_ERROR << "Dimension mismatch in gemm : A : "<< Ca().shape() <<" while B : "<<Cb().shape();
    char trans_a= get_trans(Ca(), true); 
    char trans_b= get_trans(Cb(), true); 
    int m1 = get_n_rows(Ca()), m2 = get_n_cols(Ca()), m3 = get_n_cols(Cb()); 
    int lda1 = get_lda(Ca()), lda2 = get_lda(Cb());
    //std::cerr  << " m1 .. "<< m1 << " "<< m2 << "  "<< m3 << " " << lda1 << " "<< lda2 << trans_a<<trans_b<< std::endl ;
    f77::gemm(&trans_a,&trans_b,m1,m3,m2,alpha, Ca().data_start(), lda1, Cb().data_start(), lda2, beta, Cc().data_start(), m1);
   }
   else { 
    const_qcache<MatrixType1> Ca(A);
    const_qcache<MatrixType2> Cb(B);
    if (!(Ca().dim1() == Cb().dim0())) TRIQS_RUNTIME_ERROR << "Dimension mismatch in gemm : A : "<< Ca().shape() <<" while B : "<<Cb().shape();
    char trans_a= get_trans(Ca(), false); 
    char trans_b= get_trans(Cb(), false); 
    int m1 = get_n_rows(Ca()), m2 = get_n_cols(Ca()), m3 = get_n_cols(Cb()); 
    int lda1 = get_lda(Ca()), lda2 = get_lda(Cb());
    f77::gemm(&trans_a,&trans_b,m1,m3,m2,alpha, Ca().data_start(), lda1, Cb().data_start(), lda2, beta, Cc().data_start(), m1);
   }
  }

 /*
 // generic version for non lapack 
 template<typename MatrixType1, typename MatrixType2, typename MatrixTypeOut> 
 typename std::enable_if< !(is_blas_lapack_type<typename MatrixType1::value_type>::value && have_same_value_type< MatrixType1, MatrixType2, MatrixTypeOut>::value) >::type 
 gemm (typename MatrixType1::value_type alpha, MatrixType1 const & A, MatrixType2 const & B, typename MatrixType1::value_type beta, MatrixTypeOut & C) { 
 gemm_generic(alpha,A,B,beta,C);
 }

 // make the generic version for non lapack 
 template<typename MatrixType, typename MatrixTypeOut> 
 void gemm_generic  (typename MatrixType::value_type alpha, MatrixType const & A, MatrixType const & B, typename MatrixType::value_type beta, MatrixTypeOut & C) { 

 // first resize if necessary and possible 
 resize_or_check_if_view(C,make_shape(A.dim0(),B.dim1()));

 // ecrire a la main 
 }

*/
}}}// namespace


#endif

