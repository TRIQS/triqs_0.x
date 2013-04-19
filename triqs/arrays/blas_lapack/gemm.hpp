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
#include "./tools.hpp"
#include "./qcache.hpp"

namespace triqs { namespace arrays { namespace blas { 
 
 using namespace blas_lapack_tools;
 namespace f77 { // overload

  extern "C" { 
   void TRIQS_FORTRAN_MANGLING(dgemm) (char *, char *, const int & , const int & , const int & , const double &, 
     const double[], const int &, const double[], const int &, const double &, double[], const int & ); 
   void TRIQS_FORTRAN_MANGLING(zgemm) (char *, char *, const int & , const int & , const int & , const std::complex<double> &,
     const std::complex<double>[], const int &, const std::complex<double>[], const int &, const std::complex<double> &, std::complex<double>[], const int & );
  }

  inline void gemm (char trans_a, char trans_b, const int & M, const int & N, const int & K, const double & alpha, 
    const double* A, const int & LDA, const double* B, const int & LDB, const double & beta, double* C, const int & LDC) { 
   TRIQS_FORTRAN_MANGLING(dgemm)(&trans_a,&trans_b,M,N,K,alpha, A, LDA, B, LDB, beta, C, LDC);
  }

  typedef std::complex<double> dcomplex;
  inline void gemm (char trans_a, char trans_b, const int & M, const int & N, const int & K, const dcomplex & alpha, 
    const dcomplex* A, const int & LDA, const dcomplex* B, const int & LDB, const dcomplex & beta, dcomplex* C, const int &  LDC) { 
   TRIQS_FORTRAN_MANGLING(zgemm)(&trans_a,&trans_b,M,N,K,alpha, A, LDA, B, LDB, beta, C, LDC);
  }
 }

 /**
  * Calls gemm on a matrix or view
  * Takes care of making temporary copies if necessary
  */
 template<typename MT1, typename MT2, typename MTOut> 
  typename std::enable_if< is_blas_lapack_type<typename MT1::value_type>::value && have_same_value_type< MT1, MT2, MTOut>::value >::type 
  gemm (typename MT1::value_type alpha, MT1 const & A, MT2 const & B, typename MT1::value_type beta, MTOut & C) { 

   // first resize if necessary and possible 
   resize_or_check_if_view(C,make_shape(A.dim0(),B.dim1()));

   // now we use qcache instead of the matrix to make a copy if necessary ... 
   // not optimal : if stride == 1, N ---> use LDA parameters
   // change the condition in the qcache construction....
   reflexive_qcache<MTOut> Cc(C);

   if (C.memory_layout_is_c()) {
    // then tC = tB tA ! 
    const_qcache<MT1> Cb(A); // note the inversion  A <-> B
    const_qcache<MT2> Ca(B); // note the inversion  A <-> B
    if (!(Ca().dim0() == Cb().dim1())) TRIQS_RUNTIME_ERROR << "Dimension mismatch in gemm : A : "<< Ca().shape() <<" while B : "<<Cb().shape();
    char trans_a= get_trans(Ca(), true); 
    char trans_b= get_trans(Cb(), true); 
    int m = (trans_a == 'N' ? get_n_rows(Ca()) : get_n_cols(Ca()));
    int n = (trans_b == 'N' ? get_n_cols(Cb()) : get_n_rows(Cb()));
    int k = (trans_a == 'N' ? get_n_cols(Ca()) : get_n_rows(Ca()));
    //std::cerr<< " about to call GEMM"<< std::endl ;
    //std::cerr<< "A = "<< Ca().shape()<< Ca()<< std::endl;
    //std::cerr<< "B = "<< Cb().shape()<< Cb()<< std::endl;
    //std::cerr<< "C c" << Cc().shape() << Cc().indexmap().strides() << std::endl;
    //std::cerr<<Ca().memory_layout_is_c() <<Ca().memory_layout_is_fortran()<<std::endl;
    //std::cerr<< get_n_rows(Ca())<<get_n_cols(Cb())<<get_n_cols(Ca()) << std::endl ;
    f77::gemm(trans_a,trans_b,m,n,k,
      alpha, Ca().data_start(), get_ld(Ca()) , Cb().data_start(), get_ld(Cb()), beta, Cc().data_start(), get_ld(Cc()));
    //std::cerr << " gemm ok "<< std::endl ; 
   }
   else {
    const_qcache<MT1> Ca(A);
    const_qcache<MT2> Cb(B);
    if (!(Ca().dim1() == Cb().dim0())) TRIQS_RUNTIME_ERROR << "Dimension mismatch in gemm : A : "<< Ca().shape() <<" while B : "<<Cb().shape();
    char trans_a= get_trans(Ca(), false); 
    char trans_b= get_trans(Cb(), false); 
    int m = (trans_a == 'N' ? get_n_rows(Ca()) : get_n_cols(Ca()));
    int n = (trans_b == 'N' ? get_n_cols(Cb()) : get_n_rows(Cb()));
    int k = (trans_a == 'N' ? get_n_cols(Ca()) : get_n_rows(Ca()));
    f77::gemm(trans_a,trans_b,m,n,k,
      alpha, Ca().data_start(), get_ld(Ca()) , Cb().data_start(), get_ld(Cb()), beta, Cc().data_start(), get_ld(Cc()));
   }
  }

 // make the generic version for non lapack types or more complex types
 // largely suboptimal 
 template<typename MT1, typename MT2, typename MTOut> 
  void gemm_generic  (typename MT1::value_type alpha, MT1 const & A, MT2 const & B, typename MT1::value_type beta, MTOut & C) { 
   // first resize if necessary and possible 
   resize_or_check_if_view(C,make_shape(A.dim0(),B.dim1()));
   if (A.dim1() != B.dim0()) TRIQS_RUNTIME_ERROR << "gemm generic : dimension mismatch "<< A.shape() << B.shape();
   C() = 0;
   for (int i=0; i<A.dim0(); ++i)
    for (int k=0; k<A.dim1(); ++k)
     for (int j=0; j<B.dim1(); ++j)
      C(i,j) += A(i,k)*B(k,j);
  }

 // generic version for non lapack 
 template<typename MT1, typename MT2, typename MTOut> 
  typename std::enable_if< !(is_blas_lapack_type<typename MT1::value_type>::value && have_same_value_type< MT1, MT2, MTOut>::value) >::type 
  gemm (typename MT1::value_type alpha, MT1 const & A, MT2 const & B, typename MT1::value_type beta, MTOut & C) { 
   gemm_generic(alpha,A,B,beta,C);
  }

}}}// namespace


#endif

